"""
Tools for building a reference set
"""
import collections
import csv
import itertools
import logging
import operator
import sqlite3
import tempfile

from romperroom import uclust
from Bio import SeqIO

from . import wrap

_ntf = tempfile.NamedTemporaryFile

# Utility stuff
def dedup_info_to_counts(fp):
    """
    Convert a guppy dedup file (seqid1, seqid2, count) into a dictionary
    mapping from seqid1->total count
    """
    result = collections.defaultdict(float)
    rows = csv.reader(fp)
    for i, _, c in rows:
        result[i] += float(c)
    return result

def _load_cluster_info(fp, header=True):
    r = csv.reader(fp)
    if header:
        # Skip
        next(r)
    return {k:v for k, v in r}

class SingletonDefaultDict(dict):
    """
    Dictionary-like object that returns the same value, regardless of key
    """
    def __init__(self, val=None):
        self.val = val

    def __getitem__(self, key):
        return self.val

_PARAMS = dict([('fasta_file', str),
       ('ref_fasta', str),
       ('ref_meta', str),
       ('ref_cluster_names', str),
       ('search_id', float),
       ('maxaccepts', int),
       ('maxrejects', int)])

def _load_params(con):
    cursor = con.cursor()
    cursor.execute('select key, val from params')
    result = {}
    for k, v in cursor:
        if v:
            v = _PARAMS[k](v)
        result[k] = v
    return result

def _table_exists(con, table_name):
    cursor = con.cursor()
    cursor.execute("""SELECT tbl_name FROM sqlite_master
WHERE type = 'table' AND tbl_name = ?""", [table_name])
    return cursor.fetchone() is not None

def _search(con, quiet=True):
    """
    Search the sequences in a file against a reference database
    """
    p = _load_params(con)

    cursor = con.cursor()
    count = 0
    ref_name = p['ref_fasta']
    with open(p['ref_cluster_names']) as fp:
        cluster_info = _load_cluster_info(fp)

    hit_cache = {}
    def add_hit(hit_name):
        ins = "INSERT INTO ref_seqs(name, cluster_name) VALUES (?, ?)"
        try:
            return hit_cache[hit_name]
        except KeyError:
            cluster = cluster_info[hit_name]
            cursor.execute(ins, [hit_name, cluster])
            hit_cache[hit_name] = cursor.lastrowid
            return cursor.lastrowid

    seq_id_cache = {}
    def get_seq_id(name):
        try:
            return seq_id_cache[name]
        except KeyError:
            cursor.execute('SELECT sequence_id FROM sequences WHERE name = ?', [name])
            v = cursor.fetchone()[0]
            seq_id_cache[name] = v
            return v

    with _ntf(prefix='usearch') as uc_fp:
        uclust.search(ref_name, p['fasta_file'], uc_fp.name, pct_id=0.9,
                trunclabels=True, maxaccepts=p['maxaccepts'],
                maxrejects=p['maxrejects'], quiet=quiet)

        records = uclust.parse_uclust_out(uc_fp)
        records = (i for i in records
                   if i.type == 'H' and i.pct_id >= p['search_id'] * 100.0)
        by_seq = uclust.hits_by_sequence(records)

        sql = """
INSERT INTO best_hits (sequence_id, hit_idx, ref_id, pct_id)
VALUES (?, ?, ?, ?)
"""
        for _, hits in by_seq:
            for i, h in enumerate(hits):
                # Hit id
                hit_id = add_hit(h.target_label)
                cursor.execute(sql, [get_seq_id(h.query_label), i, hit_id, h.pct_id])
                count += 1

    return count

def _load_sequences(con, sequence_file, weights=None):
    """
    Load sequences from sequence_file into database
    """
    if weights is None:
        weights = SingletonDefaultDict(1.0)

    sql = """INSERT INTO sequences (name, length, weight)
VALUES (?, ?, ?)"""

    sequences = SeqIO.parse(sequence_file, 'fasta')
    records = ((i.id, len(i), weights[i.id])
               for i in sequences)
    cursor = con.cursor()
    cursor.executemany(sql, records)
    return cursor.rowcount

def _create_tables(con, ref_fasta, ref_meta, ref_cluster_names, fasta_file,
        maxaccepts=1, maxrejects=8, search_id=0.99, quiet=True):
    cursor = con.cursor()
    cursor.executescript(SCHEMA)
    # Save parameters
    rows = [(k, locals().get(k)) for k in _PARAMS.keys()]
    cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)

Cluster = collections.namedtuple('Cluster', ['cluster_id', 'count', 'weight'])

class SearchedSequences(object):
    def __init__(self, con):
        self.con = con
        if not _table_exists(con, 'params'):
            raise ValueError("Missing table: 'params'")
        self.params = _load_params(con)

    def total_weight(self):
        cursor = self.con.cursor()
        return cursor.execute("SELECT SUM(weight) FROM sequences").fetchone()[0]

    def get_cluster(self, cluster_id):
        """
        Returns a Cluster object
        """
        cursor = self.con.cursor()
        cursor.execute("""SELECT cluster_id, cluster_count, total_weight
FROM merged_clusters WHERE cluster_id = ?""", [cluster_id])
        result = cursor.fetchone()
        if not result:
            raise KeyError(cluster_id)
        return Cluster(*result)

    def hits_by_cluster(self, hits_per_cluster=30, cluster_factor=1,
            max_per_seq=1000):
        """
        Generates an iterable of hits per cluster.

        Each item consists of:
        (cluster_number, cluster_prop, cluster_seqs, cluster_hits)
        """
        cursor = self.con.cursor()
        cursor.execute("""SELECT cluster_id, cluster_count FROM merged_clusters""")
        for cluster_id, count in cursor:
            select_count = hits_per_cluster + cluster_factor * (count - 1)
            # Find sequences
            yield (cluster_id, self._sequences_in_cluster(cluster_id),
                self._hits_in_cluster(cluster_id,
                    select_count, max_per_seq))

    def _sequences_in_cluster(self, cluster_id):
        cursor = self.con.cursor()
        cursor.execute("""SELECT s.sequence_id, s.name, s.weight, s.residues
FROM sequences s
WHERE s.cluster_id = ?""", [cluster_id])
        return [dict(zip(i.keys(), i)) for i in cursor]

    def _hits_in_cluster(self, cluster_id, hits_per_cluster=30, max_per_seq=1000):
        cursor = self.con.cursor()
        hit_sql = """
SELECT bh.name as best_hit_name, bh.ref_id
FROM best_hits bh
INNER JOIN sequences s USING(sequence_id)
WHERE s.cluster_id = ? AND bh.hit_idx < ?
ORDER BY bh.ref_id, bh.hit_idx, s.weight DESC, s.length DESC"""
        cursor.execute(hit_sql, [cluster_id, max_per_seq])
        r = (dict(zip(i.keys(), i)) for i in cursor)

        # Choose the top `hits_per_cluster` distinct hits.
        r = wrap.unique(r, operator.itemgetter('best_hit_name'))
        l = list(itertools.islice(r, 0, hits_per_cluster))
        return l

# Database schema
SCHEMA = """
CREATE TABLE sequences (
  sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR,
  length INT,
  weight FLOAT
);

CREATE UNIQUE INDEX ix_sequences_name ON sequences(name);

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
"""


def create_database(con, fasta_file, ref_fasta, ref_meta, ref_cluster_info, weights=None,
        maxaccepts=1, maxrejects=8, search_id=0.99,
        quiet=True):
    """
    Create a database of sequences searched against a sequence database for
    reference set creation.

    con: Database connection
    fasta_file: query sequences
    """
    con.row_factory = sqlite3.Row

    if _table_exists(con, 'params'):
        raise ValueError("Database exists")
    logging.info("Creating database")

    with con:
        _create_tables(con, maxaccepts=maxaccepts, maxrejects=maxrejects,
                search_id=search_id, quiet=quiet, ref_fasta=ref_fasta,
                ref_meta=ref_meta, ref_cluster_names=ref_cluster_info,
                fasta_file=fasta_file)

        seq_count = _load_sequences(con, fasta_file, weights=weights)
        logging.info("Inserted %d sequences", seq_count)

    with con:
        logging.info("Searching")
        _search(con, quiet=quiet)

def open_database(con):
    """
    Open con as a SearchSequences object
    """
    con.row_factory = sqlite3.Row
    return SearchedSequences(con)
