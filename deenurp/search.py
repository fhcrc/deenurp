"""
Tools for building a reference set
"""
import collections
import contextlib
import csv
import os
import os.path
import shutil
import sqlite3
import tempfile

from romperroom import uclust
from romperroom.RefsetInternalFasta import line_to_header
from Bio import SeqIO

_ntf = tempfile.NamedTemporaryFile

# Utility stuff
@contextlib.contextmanager
def temp_copy(path, **kwargs):
    with open(path) as fp, _ntf(delete=False, **kwargs) as tf:
        try:
            shutil.copyfileobj(fp, tf)
            tf.close()
            yield tf.name
        finally:
            os.remove(tf.name)

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

class DeNovoRefset(object):
    """
    Database housing metadata for reference set construction, namely:

    * Sequences used to search for a reference set
    * Best hits against a given reference database
    * Clustering information
    """
    _params = (('fasta_file', str),
               ('sequence_database', str),
               ('cluster_id', float),
               ('search_id', float),
               ('maxaccepts', int),
               ('maxrejects', int))

    def __init__(self, con, fasta_file=None, sequence_database=None, weights=None,
            maxaccepts=1, maxrejects=8, cluster_id=0.99, search_id=0.99,
            quiet=True):
        self.con = con
        self.con.row_factory = sqlite3.Row
        if not self._table_exists('params'):
            self.fasta_file = os.path.abspath(fasta_file)
            self.sequence_database = sequence_database
            self.weights = weights if weights is not None else SingletonDefaultDict(1.0)
            self.cluster_id = cluster_id
            self.search_id = search_id
            self.quiet = quiet
            self.maxaccepts = maxaccepts
            self.maxrejects = maxrejects

            self._create_tables()
        else:
            self._load_params()

    @classmethod
    def create(cls, con, sequence_file, sequence_database, **kwargs):
        i = cls(con, sequence_file, sequence_database, **kwargs)
        i._load_sequences()
        i._cluster()
        i._search_all()
        return i

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
INNER JOIN cluster_sequences cs ON cs.sequence_name = s.name
WHERE cs.cluster_id = ?""", [cluster_id])
        return [dict(zip(i.keys(), i)) for i in cursor]

    def _hits_in_cluster(self, cluster_id, hits_per_cluster=30, max_per_seq=1000):
        cursor = self.con.cursor()
        hit_sql = """
SELECT bh.name as best_hit_name, bh.tax_id, COUNT(*) AS hit_count
FROM best_hits bh
INNER JOIN sequences s USING(sequence_id)
INNER JOIN cluster_sequences cs ON cs.sequence_name = s.name
WHERE cs.cluster_id = ? AND bh.hit_idx < ?
GROUP BY bh.name, bh.tax_id
ORDER BY bh.hit_idx, s.weight DESC, s.length DESC
LIMIT ?"""
        cursor.execute(hit_sql, [cluster_id, max_per_seq, hits_per_cluster])
        return [dict(zip(i.keys(), i)) for i in cursor]

    def _save_params(self):
        cursor = self.con.cursor()
        rows = ((k, getattr(self, k)) for k, _ in self._params)
        cursor.executemany("INSERT INTO params VALUES (?, ?)", rows)

    def _load_params(self):
        cursor = self.con.cursor()
        cursor.execute('select key, val from params')
        d = dict(self._params)
        for k, v in cursor:
            v = d[k](v)
            setattr(self, k, v)

    def _create_tables(self):
        cursor = self.con.cursor()
        cursor.executescript(SCHEMA)
        self._save_params()

    def _table_exists(self, table_name):
        """
        Returns a boolean indicating whether table_name exists in the database.
        """
        cursor = self.con.cursor()
        cursor.execute("""SELECT tbl_name FROM sqlite_master
WHERE type = 'table' AND tbl_name = ?""", [table_name])
        return cursor.fetchone() is not None

    def _load_sequences(self):
        sql = """INSERT INTO sequences (name, residues, length, weight)
VALUES (?, ?, ?, ?)"""

        sequences = SeqIO.parse(self.fasta_file, 'fasta')
        records = ((i.id, str(i.seq), len(i), self.weights[i.id])
                   for i in sequences)
        cursor = self.con.cursor()
        with self.con:
            cursor.executemany(sql, records)
        return cursor.rowcount

    def _cluster(self):
        """
        Sort and cluster the sequences from fasta_file, loading results into
        the database.
        """
        fasta_file = self.fasta_file
        cursor = self.con.cursor()
        with tempfile.NamedTemporaryFile() as ntf:
            uclust.sort_and_cluster(fasta_file, ntf.name,
                    pct_id=self.cluster_id, trunclabels=True, quiet=self.quiet)
            records = uclust.parse_uclust_out(ntf)
            records = (i for i in records if i.type in ('H', 'S'))
            with self.con:
                for record in records:
                    # Add a cluster number
                    if record.type == 'S':
                        sql = """INSERT INTO clusters (cluster_id) VALUES (?)"""
                        cursor.execute(sql, [record.cluster_number])
                    sql = """
INSERT INTO cluster_sequences (cluster_id, orig_cluster_id, sequence_name, is_seed)
VALUES (?, ?, ?, ?)"""
                    cursor.execute(sql, [record.cluster_number, record.cluster_number,
                        record.query_label, record.type == 'S'])

    def _search_all(self):
        """
        Search all sequences against sequence_database, loading the results
        into best_hits.
        """
        def parse_hit(i, h):
            """
            Parse the target_label as a RefsetInternalFasta, returning a row
            for insertion
            """
            hit_id, annotations = line_to_header(h.target_label)
            species = None
            for rank, taxid in annotations.get('lineage', []):
                if rank == 'species':
                    species = taxid
            is_type = annotations.get('is_type') == 'type'
            query = h.query_label.split(None, 1)[0]
            return (query, i, hit_id, h.pct_id, is_type, species)

        cursor = self.con.cursor()

        # USEARCH everything
        with _ntf(prefix='seqs', suffix='.fasta') as seq_fp, \
             _ntf(prefix='usearch') as uc_fp:
            self._to_fasta(seq_fp)
            seq_fp.flush()
            uclust.search(self.sequence_database, seq_fp.name, uc_fp.name,
                    pct_id=self.search_id, trunclabels=False,
                    maxaccepts=self.maxaccepts, maxrejects=self.maxrejects,
                    quiet=self.quiet)

            by_seq = uclust.hits_by_sequence(uclust.parse_uclust_out(uc_fp))

            sql = """
INSERT INTO best_hits (sequence_id, hit_idx, name, pct_id, is_type, tax_id)
VALUES (?, ?, ?, ?, ?, ?)
"""
            count = 0
            with self.con:
                for _, hits in by_seq:
                    records = (parse_hit(i, h) for i, h in enumerate(hits))
                    cursor.executemany(sql, records)
                    count += cursor.rowcount
            return count

    def _to_fasta(self, output_fp):
        cursor = self.con.cursor()
        cursor.execute("SELECT sequence_id, name, residues FROM sequences")
        count = 0
        for record in cursor:
            output_fp.write('>{0} {1}\n{2}\n'.format(*record))
            count += 1
        return count

# Database schema
SCHEMA = """
CREATE TABLE sequences (
  sequence_id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR,
  residues VARCHAR,
  length INT,
  weight FLOAT
);

CREATE UNIQUE INDEX IX_SEQUENCES_NAME ON sequences(name);

CREATE TABLE clusters (
  cluster_id INTEGER PRIMARY KEY AUTOINCREMENT
);

CREATE TABLE cluster_sequences (
  cluster_id INTEGER REFERENCES clusters(cluster_id) ON DELETE CASCADE,
  orig_cluster_id INTEGER REFERENCES clusters(cluster_id),
  sequence_name INTEGER REFERENCES sequences(name) ON DELETE CASCADE,
  is_seed TINYINT DEFAULT 0
);

CREATE INDEX IX_CLUSTER_SEQUENCES_CLUSTER_ID ON cluster_sequences(cluster_id);

CREATE TABLE best_hits (
  hit_id INTEGER PRIMARY KEY AUTOINCREMENT,
  sequence_id INTEGER REFERENCES sequences(sequence_id) ON DELETE CASCADE,
  hit_idx INT,
  name VARCHAR,
  pct_id FLOAT,
  is_type TINYINT DEFAULT 0,
  tax_id VARCHAR
);

CREATE TABLE params (
  key VARCHAR PRIMARY KEY,
  val VARCHAR
);
"""
