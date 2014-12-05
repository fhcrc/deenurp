import os
import unittest
import numpy

from deenurp import outliers
from deenurp.subcommands.filter_outliers import filter_sequences

"""
Here's how some of the data files used in testing this module came to be:

src=/shared/silo_researcher/Matsen_F/MatsenGrp/micro_refset/rdp/10_30/1200bp/named
dest=~/src/taxtastic/testfiles

csvgrep $src/rdp_10_30_named_1200bp.seq_info.csv \
-c description -m 'Enterococcus faecium' | head -101 > $dest/e_faecium.seq_info.csv
csvcut $dest/e_faecium.seq_info.csv -c 1 | \
seqmagick convert --include-from-file - $src/rdp_10_30_named_1200bp.fasta $dest/e_faecium.fasta
cmalign.py $dest/e_faecium.fasta $dest/e_faecium.sto
seqmagick convert $dest/e_faecium.sto $dest/e_faecium.aln.fasta
FastTree -nt -makematrix $dest/e_faecium.aln.fasta > $dest/e_faecium.distmat

csvgrep $src/rdp_10_30_named_1200bp.seq_info.csv \
-c description -m 'Enterococcus faecalis' | head -101 > $dest/e_faecalis.seq_info.csv
csvcut $dest/e_faecalis.seq_info.csv -c 1 | \
seqmagick convert --include-from-file - $src/rdp_10_30_named_1200bp.fasta $dest/e_faecalis.fasta
cmalign.py $dest/e_faecalis.fasta $dest/e_faecalis.sto
seqmagick convert $dest/e_faecalis.sto $dest/e_faecalis.aln.fasta
FastTree -nt -makematrix $dest/e_faecalis.aln.fasta > $dest/e_faecalis.distmat
"""


def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)


class TestReadDists(unittest.TestCase):

    def test01(self):
        with open(data_path('e_faecium.distmat')) as f:
            taxa, mat = outliers.read_dists(f)
            self.assertEqual(len(taxa), 100)
            self.assertEqual(mat.shape, (100, 100))
            self.assertAlmostEqual(list(mat[(0, 0, 99, 99), (0, 99, 0, 99)]),
                                   [0, 0.003299, 0.003299, 0])


class TestFastTreeDists(unittest.TestCase):

    def test01(self):
        fa = data_path('e_faecium.aln.fasta')
        taxa, mat = outliers.fasttree_dists(fa)
        self.assertEqual(len(taxa), 100)
        self.assertEqual(mat.shape, (100, 100))


class TestFindOutliers(unittest.TestCase):

    def test01(self):
        with open(data_path('e_faecium.distmat')) as f:
            taxa, mat = outliers.read_dists(f)
            _, _, is_outlier = outliers.outliers(mat, cutoff=0.015)
            out = {t for t, o in zip(taxa, is_outlier) if o}
            self.assertEqual(len(out), 7)

    def test02(self):
        """
        Test special handling for matrices of size (2,2): these never
        fail.
        """

        mat = numpy.matrix([0, 0.02, 0.02, 0])
        mat.shape = (2, 2)
        _, _, is_outlier = outliers.outliers(mat, cutoff=0.015)
        self.assertFalse(any(is_outlier))

        _, _, is_outlier = outliers.outliers(mat, cutoff=0.025)
        self.assertFalse(any(is_outlier))


class TestFilterSequences(unittest.TestCase):

    def test01(self):
        fa = data_path('test_db_head.fasta')
        to_prune = filter_sequences(fa, '53635', 0.015, aligner='cmalign')
        self.assertEqual(sum(to_prune['is_out']), 5)

    def test02(self):
        fa = data_path('test_db_head.fasta')
        to_prune = filter_sequences(fa, '53635', 0.015, aligner='muscle')
        self.assertEqual(sum(to_prune['is_out']), 5)

    def test03(self):
        fa = data_path('test_db_head.fasta')
        to_prune = filter_sequences(fa, '53635', 0.015, aligner='usearch')
        self.assertEqual(sum(to_prune['is_out']), 5)



def suite():
    s = unittest.TestSuite()
    classes = [TestReadDists, TestFastTreeDists, TestFindOutliers]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))

    return s
