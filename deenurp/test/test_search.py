import collections
import os.path
from io import StringIO
import unittest
from deenurp import search

class RandomDict(dict):
    def __getitem__(self, key):
        if key not in self:
            self[key] = abs(hash(key) % 100)
        return super(RandomDict, self).__getitem__(key)

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

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
        self.assertEqual({'FUM0LCO02JSVPV': {'default': 1.0},
                          'FUM0LCO01DPWLL': {'default': 3.0},
                          'GA05AQR02HZ4VX': {'default': 3.0}}, actual)

    def test_with_sample_map(self):
        m = {'FUM0LCO01DPWLL': 'A',
             'GA05AQR02H61T7': 'B',
             'GA05AQR02HZ4VX': 'B',
             'GA05AQR02IGIUG': 'C',
             'GLKT0ZE02GQ2FO': 'A',
             'FUM0LCO02JSVPV': 'A'}
        actual = search.dedup_info_to_counts(self.s, m)
        self.assertEqual(
            {'FUM0LCO02JSVPV': {'A': 1.0},
             'FUM0LCO01DPWLL': {'A': 2.0, 'B': 1.0},
             'GA05AQR02HZ4VX': {'A': 1.0, 'B': 1.0, 'C': 1.0}},
            actual)

TestHit = collections.namedtuple('TestHit', ['query_label', 'target_label', 'pct_id'])


class SelectHitsTestCase(unittest.TestCase):
    def setUp(self):
        self.items = [
            ('seq1', [TestHit('seq1', 't1', 99.8),
                      TestHit('seq1', 't2', 99.9),
                      TestHit('seq1', 't3', 99.6)]),
            ('seq2', [TestHit('seq2', 't6', 98.4)])]

    def test_basic(self):
        r = list(search.select_hits(self.items))
        expected = [
            ('seq1', [TestHit('seq1', 't2', 99.9)]),
            ('seq2', [TestHit('seq2', 't6', 98.4)])]
        self.assertItemsEqual(expected, r)

