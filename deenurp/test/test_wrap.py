import operator
import os.path
import unittest

from Bio import SeqIO

from deenurp import wrap

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

class CmAlignTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = list(SeqIO.parse(data_path('test_input.fasta'), 'fasta'))

    def test_nompi(self):
        result = list(wrap.cmalign(self.sequences))
        self.assertEqual(len(self.sequences), len(result))

    def test_mpi(self):
        result = list(wrap.cmalign(self.sequences, mpi_args=['-np', '2']))
        self.assertEqual(len(self.sequences), len(result))


@unittest.skip("Running voronoi, rather than vorotree")
class VorotreeTestCase(unittest.TestCase):
    def setUp(self):
        self.tree_path = data_path('small.tre')

    def do_test(self, algorithm):
        actual = wrap.vorotree(self.tree_path, 2, algorithm)
        self.assertEqual(2, len(actual))

    def test_greedy(self):
        self.do_test('greedy')

    def test_full(self):
        self.do_test('full')

    @unittest.skip("Not implemented in rppr")
    def test_pam(self):
        self.do_test('pam')

class UniqueTestCase(unittest.TestCase):
    def test_nokey(self):
        l = [1, 2, 3, 1, 2, 4, 6, 7, 5, 1]
        expected = [1, 2, 3, 4, 6, 7, 5]
        actual = list(wrap.unique(l))
        self.assertEqual(expected, actual)

    def test_key(self):
        keys = ('n', 's')
        v = [(1, 'test'), (2, 'test'), (2, 'other')]
        l = [dict(zip(keys, i)) for i in v]
        expected1 = [{'n': 1, 's': 'test'},
                    {'n': 2, 's': 'test'}]
        actual1 = wrap.unique(l, key=operator.itemgetter('n'))
        self.assertEqual(expected1, list(actual1))
