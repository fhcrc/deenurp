import os.path
import sqlite3
from cStringIO import StringIO
import unittest
from deenurp import search

class RandomDict(dict):
    def __getitem__(self, key):
        if key not in self:
            self[key] = abs(hash(key) % 100)
        return super(RandomDict, self).__getitem__(key)

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

class RefsetDatabaseTestCase(unittest.TestCase):
    def setUp(self):
        self.con = sqlite3.connect(':memory:')
        self.cursor = self.con.cursor()
        self.seq_path = data_path('test_input.fasta')
        self.db_path = data_path('test_db.fasta')
        self.db = search.DeNovoRefset(self.con, self.seq_path, self.db_path,
                quiet=True)

    def tearDown(self):
        self.con.close()

    def test_create(self):
        """Ensure all tables are created"""
        expected = frozenset(['sequences', 'clusters', 'cluster_sequences',
            'best_hits'])
        for i in expected:
            self.assertTrue(self.db._table_exists(i), msg=i)

    def test_load_seqs(self):
        self.assertEqual(10, self.db._load_sequences())

    def test_load_weighted_seqs(self):
        d = RandomDict()
        self.db.weights = d
        self.assertEqual(10, self.db._load_sequences())
        s, = self.cursor.execute("""SELECT SUM(weight) FROM sequences""").fetchone()
        self.assertEqual(sum(d.values()), s)

    def test_cluster(self):
        self.db._load_sequences()
        self.db._cluster()
        expected = [(0, u'S000438419', 1), (0, u'S000871964', 0), (1,
            u'S000887598', 1), (0, u'S001610627', 0), (2, u'S002287639', 1),
            (3, u'S000136473', 1), (4, u'S000137243', 1), (5, u'S001416053',
                1), (6, u'S002222525', 1), (7, u'S000750001', 1)]
        self.assertEqual(expected, list(self.cursor.execute("SELECT * FROM cluster_sequences")))

    def test_search_each_cluster(self):
        self.db._load_sequences()
        self.db._cluster()
        result = self.db._search_all()
        self.assertEqual(10, result)

class DedupInfoToCountsTestCase(unittest.TestCase):
    def setUp(self):
        self.s = StringIO("""FUM0LCO02JSVPV,FUM0LCO02JSVPV,1
FUM0LCO01DPWLL,FUM0LCO01DPWLL,2
FUM0LCO01DPWLL,GA05AQR02H61T7,1
GA05AQR02HZ4VX,GA05AQR02HZ4VX,1
GA05AQR02HZ4VX,GA05AQR02IGIUG,1
GA05AQR02HZ4VX,GLKT0ZE02GQ2FO,1""")

    def test_basic(self):
        actual = search.dedup_info_to_counts(self.s)
        self.assertEqual({'FUM0LCO02JSVPV': 1.0,
            'FUM0LCO01DPWLL': 3.0,
            'GA05AQR02HZ4VX': 3.0}, actual)
