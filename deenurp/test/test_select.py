import os.path
import unittest

from Bio import SeqIO

from deenurp import select

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

class CmAlignTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = list(SeqIO.parse(data_path('test_input.fasta'), 'fasta'))

    def test_nompi(self):
        result = list(select.cmalign(self.sequences))
        self.assertEqual(len(self.sequences), len(result))

    def test_mpi(self):
        result = list(select.cmalign(self.sequences, mpi_args=['-np', '2']))
        self.assertEqual(len(self.sequences), len(result))


class VorotreeTestCase(unittest.TestCase):
    def setUp(self):
        self.tree_path = data_path('small.tre')

    def do_test(self, algorithm):
        actual = select.vorotree(self.tree_path, 2, algorithm)
        self.assertEqual(2, len(actual))

    def test_greedy(self):
        self.do_test('greedy')

    def test_full(self):
        self.do_test('full')

    @unittest.skip("Not implemented in rppr")
    def test_pam(self):
        self.do_test('pam')

class MergeClustersTestCase(unittest.TestCase):
    def test_multi_overlap(self):
        inp = [(1, [1, 2, 3]),
               (4, [9, 10]),
               (5, [15, 16]),
               (6, [16, 18]),
               (8, [3, 18])]
        expected = frozenset([frozenset([1, 2, 3, 18, 16, 15]), frozenset([9, 10])])
        self.assertEqual(expected, select.merge_clusters(inp))

    def test_no_overlap(self):
        inp = [(1, [1, 2]),
               (2, [4, 8])]
        expected = frozenset(frozenset(j) for i, j in inp)
        self.assertEqual(expected, select.merge_clusters(inp))

    def test_single_overlap(self):
        inp = [(1, [1, 2]),
               (2, [4, 8]),
               (3, [1, 16])]
        expected = frozenset([frozenset([1, 2, 16]), frozenset([4, 8])])
        self.assertEqual(expected, select.merge_clusters(inp))
