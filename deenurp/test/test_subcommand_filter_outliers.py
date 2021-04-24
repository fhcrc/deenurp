import unittest

import pandas as pd
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

        # confirm pairwise comparisons in the file
        tab = pd.read_table(filename, header=None, names=filter_outliers.BLAST6NAMES)
        for __, row in tab.iterrows():
            dist = distmat[seqnames.index(row['query']),
                           seqnames.index(row['target'])]
            self.assertAlmostEqual(dist, 1 - (row['pct_id'] / 100.0))


    @unittest.skipUnless(which(wrap.VSEARCH), "{} not found.".format(wrap.VSEARCH))
    def test_distmat_pairwise_vsearch(self):
        infile = util.data_path('e_faecalis.head.fasta')
        taxa, distmat = filter_outliers.distmat_pairwise(
            infile, 'foo', 'vsearch', wrap.VSEARCH)
        self.assertEqual(distmat.shape[0], len(taxa))
