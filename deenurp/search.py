"""
Tools for building a reference set
"""
import collections
import csv
import logging
import operator
import sqlite3
import tempfile

from romperroom import uclust
from Bio import SeqIO

_ntf = tempfile.NamedTemporaryFile

from .util import SingletonDefaultDict, memoize

SELECT_THRESHOLD = 0.05

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

def _load_cluster_info(fp):
    r = csv.DictReader(fp)
    return {i['seqname']: i['cluster'] for i in r}

# Parameters stores in the `params` table, with types
_PARAMS = dict([('fasta_file', str),
       ('ref_fasta', str),
       ('ref_meta', str),
       ('search_id', float),
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
    cursor.execute("""SELECT tbl_name FROM sqlite_master
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

def _search(con, quiet=True, select_threshold=SELECT_THRESHOLD):
    """
    Search the sequences in a file against a reference database
    """
    p = load_params(con)

    cursor = con.cursor()
    count = 0
    ref_name = p['ref_fasta']
    with open(p['ref_meta']) as fp:
        cluster_info = _load_cluster_info(fp)

    @memoize
    def add_hit(hit_name):
        ins = "INSERT INTO ref_seqs(name, cluster_name) VALUES (?, ?)"
        cluster = cluster_info[hit_name]
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

def _create_tables(con, ref_fasta, ref_meta, fasta_file,
        maxaccepts=1, maxrejects=8, search_id=0.99, quiet=True):
    cursor = con.cursor()
    cursor.executescript(SCHEMA)
    # Save parameters
    rows = [(k, locals().get(k)) for k in _PARAMS.keys()]
    cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)

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

CREATE VIEW vw_cluster_weights AS
SELECT cluster_name, SUM(weight) AS total_weight FROM
(SELECT DISTINCT s.sequence_id, s.weight, ref_seqs.cluster_name
 FROM sequences s
     INNER JOIN best_hits USING (sequence_id)
     INNER JOIN ref_seqs USING (ref_id)) q
GROUP BY cluster_name;
"""

def create_database(con, fasta_file, ref_fasta, ref_meta, weights=None,
        maxaccepts=1, maxrejects=8, search_id=0.99, select_threshold=SELECT_THRESHOLD,
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
                ref_meta=ref_meta, fasta_file=fasta_file)

        seq_count = _load_sequences(con, fasta_file, weights=weights)
        logging.info("Inserted %d sequences", seq_count)

    with con:
        logging.info("Searching")
        _search(con, quiet=quiet, select_threshold=select_threshold)
