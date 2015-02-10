import unittest

from Bio import SeqIO

from deenurp import wrap
from deenurp.subcommands import filter_outliers
from deenurp.util import which

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

    @unittest.skipUnless(which(wrap.VSEARCH), "{} not found.".format(wrap.VSEARCH))
    def test_distmat_pairwise_vsearch(self):
        infile = util.data_path('e_faecalis.head.fasta')
        taxa, distmat = filter_outliers.distmat_pairwise(
            infile, 'foo', 'vsearch', wrap.VSEARCH)
        self.assertEqual(distmat.shape[0], len(taxa))
