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

class SingletonDefaultDict(dict):
    """
    Dictionary-like object that returns the same value, regardless of key
    """
    def __init__(self, val=None):
        self.val = val

    def __getitem__(self, key):
        return self.val

_PARAMS = dict([('fasta_file', str),
       ('sequence_database', str),
       ('cluster_id', float),
       ('search_id', float),
       ('maxaccepts', int),
       ('maxrejects', int)])

def _load_params(con):
    cursor = con.cursor()
    cursor.execute('select key, val from params')
    result = {}
    for k, v in cursor:
        v = _PARAMS[k](v)
        result[k] = v
    return result

def _table_exists(con, table_name):
    cursor = con.cursor()
    cursor.execute("""SELECT tbl_name FROM sqlite_master
WHERE type = 'table' AND tbl_name = ?""", [table_name])
    return cursor.fetchone() is not None

def _cursor_to_fasta(cursor, output_fp):
    count = 0
    for record in cursor:
        output_fp.write('>{0}  {1}\n{2}\n'.format(*record))
        count += 1
    return count

def _unsearched_to_fasta(con, output_fp):
    """
    Write all sequences from clusters with no hits to a FASTA file
    """
    cursor = con.cursor()

    cursor.execute("""SELECT sequence_id, name, residues FROM sequences
WHERE cluster_id NOT IN (SELECT DISTINCT cluster_id
                         FROM sequences
                         INNER JOIN best_hits USING (sequence_id));""")
    return _cursor_to_fasta(cursor, output_fp)


def _desc_weight_to_fasta(con, output_fp):
    """
    Sort sequences by descending weight, then length; write to output_fp
    """
    cursor = con.cursor()
    cursor.execute("""SELECT sequence_id, name, residues
FROM sequences ORDER BY weight DESC, length DESC""")
    return _cursor_to_fasta(cursor, output_fp)


def _cluster(con, sequence_file, quiet=True):
    """
    Sort and cluster the sequences from fasta_file, loading results into
    the database.
    """
    cluster_id = _load_params(con)['cluster_id']
    cursor = con.cursor()
    with tempfile.NamedTemporaryFile() as ntf, \
         tempfile.NamedTemporaryFile(suffix='.fasta') as fa:
        _desc_weight_to_fasta(con, fa)
        fa.flush()
        uclust.cluster(fa.name, ntf.name, pct_id=cluster_id,
                trunclabels=True, quiet=quiet, usersort=True)
        records = uclust.parse_uclust_out(ntf)
        records = (i for i in records if i.type in ('H', 'S'))

        # Insert
        for record in records:
            # Add a cluster number
            if record.type == 'S':
                sql = """INSERT INTO clusters (cluster_id) VALUES (?)"""
                cursor.execute(sql, [record.cluster_number])
            sql = """
            UPDATE sequences SET cluster_id = :1, orig_cluster_id = :1
            WHERE sequence_id = :2"""
            cursor.execute(sql, [record.cluster_number, record.query_label])
            assert cursor.rowcount == 1

# Merging
def _merge_by_hit(it):
    """
    Given an iterable of (best_hit, [cluster_id, [cluster_id...]]) pairs, returns
    sets of clusters to merge based on shared best_hits.
    """
    d = {}
    for best_hit, cluster_ids in it:
        s = set(cluster_ids)

        # For each cluster, add any previous clusterings to the working cluster
        for cluster_id in cluster_ids:
            if cluster_id in d:
                s |= d[cluster_id]

        s = frozenset(s)
        for cluster_id in s:
            d[cluster_id] = s

    return frozenset(d.values())

def _find_clusters_to_merge(con, max_hits=None):
    """
    Find clusters to merge

    con: Database connection
    max_hits: number of hits to consider in each cluster [default: consider all hits]
    """
    cursor = con.cursor()
    # Generate a temporary table
    cursor.execute("""
CREATE TEMPORARY TABLE hits_to_cluster
AS
SELECT DISTINCT b.name as best_hit, cluster_id
FROM best_hits b
INNER JOIN sequences s USING(sequence_id)
{0};
""".format("" if not max_hits else "WHERE hit_idx < ?"),
    [] if not max_hits else [max_hits])
    cursor.execute("CREATE INDEX IX_hits_to_cluster_best_hit ON hits_to_cluster(best_hit);")

    try:
        sql = """
SELECT cluster_id, best_hit
FROM hits_to_cluster
WHERE best_hit IN (SELECT best_hit
                   FROM hits_to_cluster
                   GROUP BY best_hit HAVING COUNT(cluster_id) > 1)
ORDER BY best_hit;
"""
        cursor.execute(sql)
        for g, v in itertools.groupby(cursor, operator.itemgetter(1)):
            yield g, [i for i, _ in v]
    finally:
        cursor.execute("""DROP TABLE hits_to_cluster""")

def _merge_clusters(con):
    """
    Merge clusters which share a common best-hit
    """
    cmd = """UPDATE sequences
SET cluster_id = ?
WHERE cluster_id = ?"""

    cursor = con.cursor()
    to_merge = _find_clusters_to_merge(con)
    merged = _merge_by_hit(to_merge)
    for merge_group in merged:
        first = min(merge_group)
        logging.info("Merging %d clusters to %d", len(merge_group), first)
        rows = ((first, i) for i in merge_group if i != first)
        cursor.executemany(cmd, rows)

def _search_all(con, sequence_databases, quiet=True):
    """
    Search all sequences against sequence_databases, loading the results
    into best_hits.

    sequence_databases are searched in order. Any sequence without a best hit
    in sequence_database[i-1] is searched in sequence_database[i].
    """
    p = _load_params(con)

    cursor = con.cursor()
    count = 0
    # USEARCH everything
    for sequence_database in sequence_databases:
        with _ntf(prefix='seqs', suffix='.fasta') as seq_fp, \
             _ntf(prefix='usearch') as uc_fp:

            # Add a reference database
            cursor.execute("""INSERT INTO refs (file_path) VALUES (?)""",
                    [sequence_database])
            ref_id = cursor.lastrowid

            u_count = _unsearched_to_fasta(con, seq_fp)
            logging.info("%d sequences to search against %s", u_count, sequence_database)
            seq_fp.flush()
            uclust.search(sequence_database, seq_fp.name, uc_fp.name,
                    pct_id=p['search_id'], trunclabels=True,
                    maxaccepts=p['maxaccepts'], maxrejects=p['maxrejects'],
                    quiet=quiet, search_pct_id=0.9)

            by_seq = uclust.hits_by_sequence(uclust.parse_uclust_out(uc_fp))

            sql = """
INSERT INTO best_hits (sequence_id, hit_idx, name, pct_id, ref_id)
VALUES (?, ?, ?, ?, ?)
"""

            for _, hits in by_seq:
                records = ((h.query_label, i, h.target_label, h.pct_id, ref_id)
                           for i, h in enumerate(hits))
                cursor.executemany(sql, records)
                count += cursor.rowcount

    return count

def _load_sequences(con, sequence_file, weights=None):
    """
    Load sequences from sequence_file into database
    """
    if weights is None:
        weights = SingletonDefaultDict(1.0)

    sql = """INSERT INTO sequences (name, residues, length, weight)
VALUES (?, ?, ?, ?)"""

    sequences = SeqIO.parse(sequence_file, 'fasta')
    records = ((i.id, str(i.seq), len(i), weights[i.id])
               for i in sequences)
    cursor = con.cursor()
    cursor.executemany(sql, records)
    return cursor.rowcount

def _create_tables(con, maxaccepts=1, maxrejects=8, cluster_id=0.99,
        search_id=0.99, quiet=True):
    cursor = con.cursor()
    cursor.executescript(SCHEMA)
    # Save parameters
    rows = [(k, locals().get(k)) for k in _PARAMS.keys()]
    cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)

class SearchedSequences(object):
    def __init__(self, con):
        self.con = con
        if not _table_exists(con, 'params'):
            raise ValueError("Missing table: 'params'")
        self.params = _load_params(con)

    def total_weight(self):
        cursor = self.con.cursor()
        return cursor.execute("SELECT SUM(weight) FROM sequences").fetchone()[0]

    def hits_by_cluster(self, hits_per_cluster=30, max_per_seq=1000):
        """
        Generates an iterable of hits per cluster.

        Each item consists of:
        (cluster_number, cluster_prop, cluster_seqs, cluster_hits)
        """
        cursor = self.con.cursor()
        cursor.execute("""SELECT cluster_id FROM clusters""")
        clusters = [i for i, in cursor]
        for cluster_id in clusters:
            # Find sequences
            yield (cluster_id, self._sequences_in_cluster(cluster_id),
                self._hits_in_cluster(cluster_id, hits_per_cluster, max_per_seq))

    def _sequences_in_cluster(self, cluster_id):
        cursor = self.con.cursor()
        cursor.execute("""SELECT s.sequence_id, s.name, s.weight, s.residues
FROM sequences s
WHERE s.cluster_id = ?""", [cluster_id])
        return [dict(zip(i.keys(), i)) for i in cursor]

    def _hits_in_cluster(self, cluster_id, hits_per_cluster=30, max_per_seq=1000):
        cursor = self.con.cursor()
        hit_sql = """
SELECT bh.name as best_hit_name, bh.ref_id, COUNT(*) AS hit_count
FROM best_hits bh
INNER JOIN sequences s USING(sequence_id)
WHERE s.cluster_id = ? AND bh.hit_idx < ?
GROUP BY bh.name, bh.ref_id
ORDER BY bh.hit_idx, s.weight DESC, s.length DESC
LIMIT ?"""
        cursor.execute(hit_sql, [cluster_id, max_per_seq, hits_per_cluster])
        return [dict(zip(i.keys(), i)) for i in cursor]

# Database schema
SCHEMA = """
CREATE TABLE clusters (
  cluster_id INTEGER PRIMARY KEY AUTOINCREMENT
);

CREATE TABLE sequences (
  sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR,
  residues VARCHAR,
  length INT,
  weight FLOAT,
  cluster_id INTEGER REFERENCES clusters(cluster_id),
  orig_cluster_id INTEGER REFERENCES clusters(cluster_id)
);

CREATE UNIQUE INDEX ix_sequences_name ON sequences(name);
CREATE INDEX ix_sequences_cluster_id ON sequences(cluster_id);

CREATE TABLE refs (
  ref_id INTEGER PRIMARY KEY AUTOINCREMENT,
  file_path VARCHAR
);

CREATE TABLE best_hits (
  hit_id INTEGER PRIMARY KEY AUTOINCREMENT,
  sequence_id INTEGER REFERENCES sequences(sequence_id) ON DELETE CASCADE,
  hit_idx INT,
  name VARCHAR,
  pct_id FLOAT,
  ref_id INT REFERENCES refs(ref_id) NOT NULL
);

CREATE INDEX ix_best_hits_sequence_id ON best_hits(sequence_id);

CREATE TABLE params (
  key VARCHAR PRIMARY KEY,
  val VARCHAR
);

CREATE VIEW cluster_hits AS
SELECT s.cluster_id, bh.name AS best_hit
FROM best_hits bh
INNER JOIN sequences s USING (sequence_id)
GROUP BY s.cluster_id, bh.name;
"""

def create_database(con, fasta_file, reference_candidates, weights=None,
        maxaccepts=1, maxrejects=8, cluster_id=0.99, search_id=0.99,
        quiet=True):
    """
    Create a database of sequences searched against a sequence database for
    reference set creation.

    con: Database connection
    fasta_file: query sequences
    reference_candidates: Iterable of file names to be searched for candidate
        reference sequences, ordered by decreasing preference. Any cluster with a
        hit in a previous database is not searched in subsequent databases.
    """
    con.row_factory = sqlite3.Row

    if _table_exists(con, 'params'):
        raise ValueError("Database exists")
    logging.info("Creating database")

    with con:
        _create_tables(con, maxaccepts=maxaccepts, maxrejects=maxrejects,
                cluster_id=cluster_id, search_id=search_id, quiet=quiet)

    with con:
        seq_count = _load_sequences(con, fasta_file, weights=weights)
        logging.info("Inserted %d sequences", seq_count)

        logging.info("Clustering")
        _cluster(con, fasta_file, quiet=quiet)

        logging.info("Searching")
        _search_all(con, reference_candidates, quiet=quiet)

    with con:
        logging.info("Merging clusters with common best-hits")
        _merge_clusters(con)

def open_database(con):
    """
    Open con as a SearchSequences object
    """
    con.row_factory = sqlite3.Row
    return SearchedSequences(con)
