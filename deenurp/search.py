"""
Tools for building a reference set
"""
import collections
import csv
import functools
import logging
import operator
import os
import sqlite3
import tempfile

from deenurp import uclust
from Bio import SeqIO

from .util import SingletonDefaultDict, memoize

_ntf = tempfile.NamedTemporaryFile

SELECT_THRESHOLD = 0.05
SEARCH_THRESHOLD = 0.90
SEARCH_IDENTITY = 0.97

# Utility stuff


def dedup_info_to_counts(fp, sample_map=None):
    """
    Convert a guppy dedup file (seqid1, seqid2, count) into a dictionary
    mapping {seqid:{sample:count}}
    """
    if sample_map is None:
        sample_map = SingletonDefaultDict('default')
    result = collections.defaultdict(
        functools.partial(
            collections.defaultdict, float))
    rows = csv.reader(fp)
    for i, j, c in rows:
        result[i][sample_map[j]] += float(c)
    return dict(result)


def load_sample_map(fp, header=False):
    """
    Load a pplacer-compatible sample map

    Returns a dict mapping from {sequence:sample}
    """
    r = csv.reader(fp)
    if header:
        next(r)
    return dict(r)


def _load_cluster_info(fp, group_field='cluster'):
    r = csv.DictReader(fp)
    return {i['seqname']: i[group_field] for i in r}

# Parameters stores in the `params` table, with types
_PARAMS = dict([('fasta_file', str),
                ('ref_fasta', str),
                ('ref_meta', str),
                ('search_identity', float),
                ('group_field', str),
                ('maxaccepts', int),
                ('maxrejects', int)])


def load_params(con):
    """
    Load parameters from the ``params`` table
    """
    cursor = con.cursor()
    sql = 'select key, val from params'
    logging.debug(sql)
    cursor.execute(sql)
    result = {}
    for k, v in cursor:
        if v:
            v = _PARAMS[k](v)
        result[k] = v
    return result


def _table_exists(con, table_name):
    """
    Returns whether or not ``table_name`` exists in ``con``
    """
    cursor = con.cursor()
    sql = ('SELECT tbl_name '
           'FROM sqlite_master '
           'WHERE type = "table" AND tbl_name = ?')
    logging.debug(sql.replace('?', '{}').format(table_name))
    cursor.execute(sql, [table_name])
    return cursor.fetchone() is not None


def select_hits(hits_by_seq, threshold=SELECT_THRESHOLD):
    """
    Select all hits for each sequence within ``threshold``
    of the best percent id
    """
    for seq, hits in hits_by_seq:
        hits = list(hits)
        hits.sort(key=operator.attrgetter('pct_id'), reverse=True)
        result = [hits[0]]
        best_pct_id = hits[0].pct_id
        result.extend(
            i for i in hits[
                1:] if best_pct_id -
            i.pct_id < threshold)
        yield seq, result


def _search(con, quiet=True, select_threshold=SELECT_THRESHOLD,
            search_threshold=SEARCH_THRESHOLD, blacklist=None):
    """
    Search the sequences in a file against a reference database
    """
    blacklist = blacklist or set()
    p = load_params(con)

    cursor = con.cursor()
    count = 0
    ref_name = p['ref_fasta']
    with open(p['ref_meta']) as fp:
        cluster_info = _load_cluster_info(fp, p['group_field'])

    @memoize
    def add_hit(hit_name, clusterj):
        ins = 'INSERT INTO ref_seqs(name, cluster_name) VALUES (?, ?)'
        logging.debug(ins.replace('?', '{}').format(hit_name, cluster))
        cursor.execute(ins, [hit_name, cluster])
        return cursor.lastrowid

    @memoize
    def get_seq_id(name):
        sql = 'SELECT sequence_id FROM sequences WHERE name = ?'
        logging.debug(sql.replace('?', '{}').format(name))
        cursor.execute(sql, [name])
        return cursor.fetchone()[0]

    with _ntf(prefix='usearch') as uc_fp:
        uclust.search(
            ref_name,
            p['fasta_file'],
            uc_fp.name,
            pct_id=search_threshold,
            maxaccepts=p['maxaccepts'],
            maxrejects=p['maxrejects'],
            quiet=quiet)

        # import shutil
        # shutil.copy(uc_fp.name, '.')

        records = uclust.parse_uclust_out(uc_fp)
        records = (i for i in records if i.type ==
                   'H' and i.pct_id >= p['search_identity'] * 100.0)
        by_seq = uclust.hits_by_sequence(records)
        by_seq = select_hits(by_seq, select_threshold)

        sql = """
INSERT INTO best_hits (sequence_id, hit_idx, ref_id, pct_id)
VALUES (?, ?, ?, ?)
"""
        for _, hits in by_seq:
            # Drop clusters from blacklist
            hits = (
                h for h in hits if not cluster_info[
                    h.target_label] in blacklist)
            seen_clusters = set()
            for i, h in enumerate(hits):
                cluster = cluster_info[h.target_label]

                # Only keep one sequence per cluster
                if cluster in seen_clusters:
                    continue
                else:
                    seen_clusters.add(cluster)

                # Hit id
                hit_id = add_hit(h.target_label, cluster)
                seq_id = get_seq_id(h.query_label)
                logging.debug(sql.replace('?', '{}').format(
                    seq_id, i, hit_id, h.pct_id))
                cursor.execute(sql, [seq_id, i, hit_id, h.pct_id])
                count += 1

    return count


def _load_sequences(con, sequence_file, weights=None):
    """
    Load sequences from sequence_file into database
    """
    if weights is None:
        weights = SingletonDefaultDict({'default': 1.0})
    seq_count = 0

    @memoize
    def get_sample_id(sample_name):
        cursor = con.cursor()
        sql = 'SELECT sample_id FROM samples WHERE name = ?'
        logging.debug(sql.replace('?', '{}').format(sample_name))
        cursor.execute(sql, [sample_name])
        result = cursor.fetchone()
        if result:
            return result[0]
        else:
            sql = """INSERT INTO samples (name) VALUES (?)"""
            logging.debug(sql.replace('?', '{}').format(sample_name))
            cursor.execute(sql, [sample_name])
            return cursor.lastrowid

    sequences = SeqIO.parse(sequence_file, 'fasta')
    cursor = con.cursor()
    sequence_insert_sql = """INSERT INTO sequences (name, length)
VALUES (?, ?)"""
    debug = sequence_insert_sql.replace('?', '{}')
    for sequence in sequences:
        seq_len = len(sequence)
        logging.debug(debug.format(sequence.id, seq_len))
        cursor.execute(sequence_insert_sql, [sequence.id, seq_len])
        seq_id = cursor.lastrowid
        seq_count += 1
        if sequence.id not in weights:
            continue
        for sample, weight in weights[sequence.id].items():
            sample_id = get_sample_id(sample)
            cursor.execute("""INSERT INTO sequences_samples
                           (sequence_id, sample_id, weight)
                           VALUES (?, ?, ?)""",
                           [seq_id, sample_id, weight])
    return seq_count


def _create_tables(
        con,
        ref_fasta,
        ref_meta,
        fasta_file,
        maxaccepts=1,
        maxrejects=8,
        search_identity=SEARCH_IDENTITY,
        quiet=True,
        group_field='cluster'):
    schema = os.path.join(os.path.dirname(__file__), 'data', 'search.schema')
    cursor = con.cursor()
    cursor.executescript(open(schema).read().strip())
    # Save parameters
    rows = [(k, locals().get(k)) for k in _PARAMS.keys()]
    cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)


def create_database(
        con,
        fasta_file,
        ref_fasta,
        ref_meta,
        weights=None,
        maxaccepts=1,
        maxrejects=8,
        search_identity=SEARCH_IDENTITY,
        select_threshold=SELECT_THRESHOLD,
        search_threshold=SEARCH_THRESHOLD,
        quiet=True,
        group_field='cluster',
        blacklist=None):
    """
    Create a database of sequences searched against a sequence database for
    reference set creation.

    con: Database connection
    fasta_file: query sequences
    search_identity:
    select_threshold:
    search_threshold:
    """

    # print "search_identity", search_identity
    # print "select_threshold", select_threshold
    # print "search_threshold", search_threshold

    # try to avoid issues with floating point comparison
    allowed_difference = 0.0001
    if search_threshold - search_identity > allowed_difference:
        msg = ('search_identity ({}) should not be less '
               'than than search_threshold ({})')
        raise ValueError(msg.format(search_identity, search_threshold))

    con.row_factory = sqlite3.Row

    blacklist = blacklist or set()

    if _table_exists(con, 'params'):
        raise ValueError("Database exists")
    logging.info("Creating database")

    with con:
        _create_tables(
            con,
            maxaccepts=maxaccepts,
            maxrejects=maxrejects,
            search_identity=search_identity,
            quiet=quiet,
            ref_fasta=ref_fasta,
            ref_meta=ref_meta,
            fasta_file=fasta_file,
            group_field=group_field)

        seq_count = _load_sequences(con, fasta_file, weights=weights)
        logging.info("Inserted %d sequences", seq_count)

    with con:
        logging.info("Searching")
        _search(con, quiet=quiet, select_threshold=select_threshold,
                search_threshold=search_threshold, blacklist=blacklist)
