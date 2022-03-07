import os
import unittest

from deenurp import outliers, wrap
from deenurp.subcommands.filter_outliers import filter_sequences, distmat_muscle
from deenurp.util import MissingDependencyError, which

try:
    import numpy as np
    import pandas as pd
except ImportError as err:
    # prefer errors within tests over failure at the time the test
    # suites are assembled
    print(err)

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


@unittest.skipUnless(which('FastTree'), "FastTree not found.")
class TestFastTreeDists(unittest.TestCase):

    def test01(self):
        fa = data_path('e_faecium.aln.fasta')
        taxa, mat = outliers.fasttree_dists(fa)
        self.assertEqual(len(taxa), 100)
        self.assertEqual(mat.shape, (100, 100))


class TestFindOutliers(unittest.TestCase):

    def setUp(self):
        with open(data_path('e_faecium.distmat')) as f:
            self.taxa, self.mat = outliers.read_dists(f)

    def test_find_medoid_01(self):
        medoid = outliers.find_medoid(self.mat)
        self.assertEqual(medoid, 3)

    def test_find_medoid_02(self):
        n, m = self.mat.shape
        medoid = outliers.find_medoid(
            self.mat, ii=np.array([i > n/2. for i in range(n)], dtype=bool))
        self.assertEqual(medoid, 82)

    def test_scipy_cluster(self):
        n, m = self.mat.shape
        clusters, title = outliers.scipy_cluster(self.mat, 'single', t=0.01)
        self.assertEqual(len(clusters), n)

    def test_find_cluster_medoids(self):
        clusters, title = outliers.scipy_cluster(self.mat, 'single', t=0.01)
        df = outliers.find_cluster_medoids(self.mat, clusters)
        self.assertEqual(len(set(clusters)), df.shape[0])

    def test_choose_clusters(self):

        s = """{"cluster":{"0":0,"1":-1,"2":1},
        "count":{"0":278,"1":17,"2":6},
        "medoid":{"0":238.0,"1":null,"2":284.0},
        "dist":{"0":0.0,"1":null,"2":0.089}}"""

        df = pd.read_json(s)
        # output is a set of cluster names (not indices)
        self.assertSetEqual(set(outliers.choose_clusters(df, 2, 0.015)), {0})
        self.assertSetEqual(set(outliers.choose_clusters(df, 2, 0.1)), {0, 1})

    def test_scaled_radius(self):
        R = outliers.scaled_radius(self.mat, percentile=90, min_radius=0.01)
        self.assertGreater(R, 0.01)
        R = outliers.scaled_radius(self.mat, percentile=90, min_radius=0.015)
        self.assertEqual(R, 0.015)
        R = outliers.scaled_radius(self.mat, percentile=90, max_radius=0.01)
        self.assertEqual(R, 0.01)

    def test_outliers_01(self):
        _, _, is_outlier = outliers.outliers(self.mat, radius=0.015)
        out = {t for t, o in zip(self.taxa, is_outlier) if o}
        self.assertEqual(len(out), 7)

    def test_outliers_02(self):
        _, _, is_outlier, _ = outliers.outliers_by_cluster(self.mat, t=0.015, D=0.015)
        out = {t for t, o in zip(self.taxa, is_outlier) if o}
        self.assertEqual(len(out), 4)

    def test_mds_01(self):
        df = outliers.mds(self.mat, self.taxa)
        self.assertEqual(df.shape[0], self.mat.shape[0])
        self.assertTrue((df['seqname'] == self.taxa).all())

    def test_mds_02(self):
        # handle all zero distances
        df = outliers.mds(np.zeros(shape=self.mat.shape), self.taxa)
        self.assertEqual(df.shape[0], self.mat.shape[0])
        self.assertTrue((df['seqname'] == self.taxa).all())
        self.assertTrue((df['x'] == 0).all())
        self.assertTrue((df['y'] == 0).all())

try:
    wrap.require_executable(wrap.VSEARCH)
except MissingDependencyError as e:
    vsearch_available = False
else:
    vsearch_available = True


@unittest.skipUnless(which('muscle') and vsearch_available,
                     "vsearch and muscle not found.")
class TestFilterSequences(unittest.TestCase):

    fa = data_path('test_db_head.fasta')
    tax_id = '53635'

    @classmethod
    def setUpClass(cls):
        cls.taxa, cls.distmat = distmat_muscle(cls.fa, '{}_'.format(cls.tax_id))

    def test01(self):
        to_prune = filter_sequences(self.tax_id,
                                    sequence_file=self.fa,
                                    strategy='radius',
                                    cutoff=0.015,
                                    aligner='cmalign')
        self.assertEqual(sum(to_prune['is_out']), 5)

    @unittest.skipUnless(vsearch_available, "{} not found.".format(wrap.VSEARCH))
    def test02(self):
        to_prune = filter_sequences(self.tax_id,
                                    sequence_file=self.fa,
                                    strategy='radius',
                                    cutoff=0.015,
                                    aligner='vsearch')
        self.assertEqual(sum(to_prune['is_out']), 5)

    def test03(self):
        to_prune = filter_sequences(self.tax_id,
                                    sequence_file=self.fa,
                                    strategy='radius',
                                    cutoff=0.015,
                                    aligner='muscle')
        self.assertEqual(sum(to_prune['is_out']), 5)

    def test04(self):
        to_prune = filter_sequences(self.tax_id,
                                    distmat=self.distmat,
                                    taxa=self.taxa,
                                    strategy='cluster',
                                    cutoff=0.015)
        self.assertEqual(sum(to_prune['is_out']), 5)

    def test05(self):
        to_prune = filter_sequences(self.tax_id,
                                    distmat=self.distmat,
                                    taxa=self.taxa,
                                    strategy='cluster',
                                    percentile=90)
        self.assertEqual(sum(to_prune['is_out']), 0)

    def test06(self):
        to_prune = filter_sequences(self.tax_id,
                                    distmat=self.distmat,
                                    taxa=self.taxa,
                                    strategy='cluster',
                                    percentile=90, min_radius=0.1)
        self.assertEqual(sum(to_prune['is_out']), 0)
