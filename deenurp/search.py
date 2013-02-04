"""
Tools for building a reference set
"""
import collections
import csv
import functools
import logging
import operator
import sqlite3
import tempfile

from deenurp import uclust
from Bio import SeqIO

_ntf = tempfile.NamedTemporaryFile

from .util import SingletonDefaultDict, memoize

SELECT_THRESHOLD = 0.05

# Utility stuff
def dedup_info_to_counts(fp, sample_map=None):
    """
    Convert a guppy dedup file (seqid1, seqid2, count) into a dictionary
    mapping {seqid:{sample:count}}
    """
    if sample_map is None:
        sample_map = SingletonDefaultDict('default')
    result = collections.defaultdict(functools.partial(collections.defaultdict, float))
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
       ('search_id', float),
       ('group_field', str),
       ('maxaccepts', int),
       ('maxrejects', int)])

def load_params(con):
    """
    Load parameters from the ``params`` table
    """
    cursor = con.cursor()
    cursor.execute('select key, val from params')
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
    cursor.execute("""SELECT tbl_name
FROM sqlite_master
WHERE type = 'table' AND tbl_name = ?""", [table_name])
    return cursor.fetchone() is not None

def select_hits(hits_by_seq, threshold=SELECT_THRESHOLD):
    """
    Select all hits for each sequence within ``threshold`` of the best percent id
    """
    for seq, hits in hits_by_seq:
        hits = list(hits)
        hits.sort(key=operator.attrgetter('pct_id'), reverse=True)
        result = [hits[0]]
        best_pct_id = hits[0].pct_id
        result.extend(i for i in hits[1:] if best_pct_id - i.pct_id < threshold)
        yield seq, result

def _search(con, quiet=True, select_threshold=SELECT_THRESHOLD, blacklist=None):
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
        ins = "INSERT INTO ref_seqs(name, cluster_name) VALUES (?, ?)"
        cursor.execute(ins, [hit_name, cluster])
        return cursor.lastrowid

    @memoize
    def get_seq_id(name):
        cursor.execute('SELECT sequence_id FROM sequences WHERE name = ?', [name])
        return cursor.fetchone()[0]

    with _ntf(prefix='usearch') as uc_fp:
        uclust.search(ref_name, p['fasta_file'], uc_fp.name, pct_id=0.9,
                trunclabels=True, maxaccepts=p['maxaccepts'],
                maxrejects=p['maxrejects'], quiet=quiet)

        records = uclust.parse_uclust_out(uc_fp)
        records = (i for i in records
                   if i.type == 'H' and i.pct_id >= p['search_id'] * 100.0)
        by_seq = uclust.hits_by_sequence(records)
        by_seq = select_hits(by_seq, select_threshold)

        sql = """
INSERT INTO best_hits (sequence_id, hit_idx, ref_id, pct_id)
VALUES (?, ?, ?, ?)
"""
        for _, hits in by_seq:
            # Drop clusters from blacklist
            hits = (h for h in hits if not cluster_info[h.target_label] in blacklist)
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
                cursor.execute(sql, [get_seq_id(h.query_label), i, hit_id, h.pct_id])
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
        cursor.execute("""SELECT sample_id FROM samples WHERE name = ?""", [sample_name])
        result = cursor.fetchone()
        if result:
            return result[0]
        else:
            cursor.execute("""INSERT INTO samples (name) VALUES (?)""", [sample_name])
            return cursor.lastrowid

    sequences = SeqIO.parse(sequence_file, 'fasta')
    cursor = con.cursor()
    sequence_insert_sql = """INSERT INTO sequences (name, length)
VALUES (?, ?)"""
    for sequence in sequences:
        cursor.execute(sequence_insert_sql, [sequence.id, len(sequence)])
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

def _create_tables(con, ref_fasta, ref_meta, fasta_file,
        maxaccepts=1, maxrejects=8, search_id=0.99, quiet=True,
        group_field='cluster'):
    cursor = con.cursor()
    cursor.executescript(SCHEMA)
    # Save parameters
    rows = [(k, locals().get(k)) for k in _PARAMS.keys()]
    cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)

# Database schema
SCHEMA = """
CREATE TABLE samples (
  sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR UNIQUE
);

CREATE TABLE sequences (
  sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR,
  length INT
);
CREATE UNIQUE INDEX ix_sequences_name ON sequences(name);

CREATE TABLE sequences_samples (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  sequence_id INT REFERENCES sequences(sequence_id) ON DELETE CASCADE,
  sample_id INT REFERENCES samples(sample_id) ON DELETE CASCADE,
  weight FLOAT
);
CREATE INDEX ix_sequences_samples_sequence_id ON sequences_samples(sequence_id);
CREATE INDEX ix_sequences_samples_sample_id ON sequences_samples(sample_id);

CREATE TABLE ref_seqs (
  ref_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR UNIQUE,
  cluster_name VARCHAR
);

CREATE INDEX ix_ref_seqs_cluster_name ON ref_seqs(cluster_name);

CREATE TABLE best_hits (
  hit_id INTEGER PRIMARY KEY AUTOINCREMENT,
  ref_id INTEGER REFERENCES ref_seqs(ref_id) ON DELETE CASCADE,
  sequence_id INTEGER REFERENCES sequences(sequence_id) ON DELETE CASCADE,
  hit_idx INT,
  pct_id FLOAT
);

CREATE INDEX ix_best_hits_sequence_id ON best_hits(sequence_id);
CREATE INDEX ix_best_hits_ref_id ON best_hits(ref_id);

CREATE TABLE params (
  key VARCHAR PRIMARY KEY,
  val VARCHAR
);

CREATE VIEW vw_cluster_weights AS
SELECT cluster_name, SUM(weight) AS total_weight FROM
(SELECT DISTINCT s.sequence_id, ss.weight as weight, ref_seqs.cluster_name
 FROM sequences s
     INNER JOIN sequences_samples ss USING (sequence_id)
     INNER JOIN best_hits USING (sequence_id)
     INNER JOIN ref_seqs USING (ref_id)) q
GROUP BY cluster_name;

CREATE VIEW vw_sample_weights AS
SELECT s.sample_id, s.name, SUM(ss.weight) AS total_weight
FROM samples s
INNER JOIN sequences_samples ss USING (sample_id)
GROUP BY s.sample_id, s.name
"""

def create_database(con, fasta_file, ref_fasta, ref_meta, weights=None,
        maxaccepts=1, maxrejects=8, search_id=0.99,
        select_threshold=SELECT_THRESHOLD, quiet=True, group_field='cluster',
        blacklist=None):
    """
    Create a database of sequences searched against a sequence database for
    reference set creation.

    con: Database connection
    fasta_file: query sequences
    """
    con.row_factory = sqlite3.Row

    blacklist = blacklist or set()

    if _table_exists(con, 'params'):
        raise ValueError("Database exists")
    logging.info("Creating database")

    with con:
        _create_tables(con, maxaccepts=maxaccepts, maxrejects=maxrejects,
                search_id=search_id, quiet=quiet, ref_fasta=ref_fasta,
                ref_meta=ref_meta, fasta_file=fasta_file, group_field=group_field)

        seq_count = _load_sequences(con, fasta_file, weights=weights)
        logging.info("Inserted %d sequences", seq_count)

    with con:
        logging.info("Searching")
        _search(con, quiet=quiet, select_threshold=select_threshold,
                blacklist=blacklist)
