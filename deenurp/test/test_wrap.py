import os.path
import subprocess
import tempfile
import unittest

from Bio import SeqIO

import deenurp
from deenurp import wrap
from deenurp.test import util
from deenurp.util import which, MissingDependencyError

@unittest.skipUnless(which('cmalign'), "cmalign not found.")
class CmAlignTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = list(SeqIO.parse(util.data_path('test_input.fasta'), 'fasta'))

    def test_oneproc(self):
        result = list(wrap.cmalign(self.sequences))
        self.assertEqual(len(self.sequences), len(result))

    def test_twoproc(self):
        result = list(wrap.cmalign(self.sequences, cpu=2))
        self.assertEqual(len(self.sequences), len(result))

    def test_defaultproc(self):
        result = list(wrap.cmalign(self.sequences))
        self.assertEqual(len(self.sequences), len(result))


class CMTestCase(unittest.TestCase):
    def test_find_cm(self):
        self.assertTrue(os.path.isfile(wrap.CM))

@unittest.skipUnless(which('FastTree'), "FastTree not found")
class AsRefpkgTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = SeqIO.parse(
            util.data_path('taxon171552.seqs.fasta'), 'fasta')

    def test_as_refpkg(self):
        with wrap.as_refpkg(self.sequences) as refpkg:
            self.assertTrue(os.path.isdir(refpkg.path))

            if which('rppr'):
                out = subprocess.check_output(['rppr', 'check', '-c', refpkg.path])
                self.assertTrue('OK!' in out, out)

@unittest.skipUnless(which('rppr'), "rppr not found")
class RpprMinAdclTreeTestCase(unittest.TestCase):
    def setUp(self):
        self.tf = tempfile.NamedTemporaryFile(prefix='adcl', suffix='.tre')
        self.tf.write("((C000721552:0.20692,C002038857:0.00015)0.844:0.01031,C002038856:0.00014,((C002963332:0.08558,(C001550734:0.06763,((C000004779:0.03889,C002963310:0.04622)0.633:0.00151,(C002963318:0.00014,C002963266:0.00014)0.697:0.00016)0.992:0.15253)0.889:0.07332)0.924:0.07668,C002038858:0.01032)0.907:0.00014);\n")
        self.tf.flush()

    def tearDown(self):
        self.tf.close()

    def test_min_adcl_prune3(self):
        prune_seqs = wrap.rppr_min_adcl_tree(self.tf.name, 7)
        self.assertEqual(3, len(prune_seqs))

    def test_min_adcl_prune7(self):
        prune_seqs = wrap.rppr_min_adcl_tree(self.tf.name, 3)
        self.assertEqual(7, len(prune_seqs))


def suite():
    s = unittest.TestSuite()
    classes = [CmAlignTestCase, CMTestCase, RpprMinAdclTreeTestCase]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))

    return s


try:
    wrap.require_executable(wrap.USEARCH)
    wrap._require_usearch6(wrap.USEARCH)
except MissingDependencyError, e:
    usearch_available = False
else:
    usearch_available = True

@unittest.skipUnless(usearch_available, "{} not found.".format(wrap.USEARCH))
class UsearchTestCase(unittest.TestCase):
    def setUp(self):
        self.sequencefile = util.data_path('test_input.fasta')

    def test_require_usearch6(self):
        wrap._require_usearch6()

    def test_usearch_allpairs_files(self):
        with deenurp.util.ntf(suffix='.blast6out') as outfile:
            wrap.usearch_allpairs_files(self.sequencefile, outfile.name)
            self.assertTrue(os.path.exists(outfile.name))

        self.assertFalse(os.path.exists(outfile.name))
