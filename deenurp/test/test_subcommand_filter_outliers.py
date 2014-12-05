import unittest

from Bio import SeqIO

from deenurp.subcommands import filter_outliers
from deenurp.test import util


class FilterOutliersFunctions(unittest.TestCase):
    def setUp(self):
        pass

    def test_parse_usearch_allpairs(self):
        filename = util.data_path('e_faecalis.head.allpairs')
        with open(util.data_path('e_faecalis.head.fasta')) as f:
            seqs = SeqIO.parse(f, 'fasta')
            seqnames = [seq.id for seq in seqs]

        distmat = filter_outliers.parse_usearch_allpairs(filename, seqnames)
        self.assertEqual(len(seqnames), distmat.shape[0])

    def test_distmat_usearch(self):
        infile = util.data_path('e_faecalis.head.fasta')
        taxa, distmat = filter_outliers.distmat_usearch(infile, 'foo')
        self.assertEqual(distmat.shape[0], len(taxa))


def suite():
    s = unittest.TestSuite()
    classes = [FilterOutliersFunctions]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))
    return s
