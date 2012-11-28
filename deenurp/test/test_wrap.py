import os.path
import subprocess
import unittest

from Bio import SeqIO

from deenurp import wrap
from deenurp.test import util

@unittest.skipUnless(util.which('cmalign'), "cmalign not found.")
class CmAlignTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = list(SeqIO.parse(util.data_path('test_input.fasta'), 'fasta'))

    def test_nompi(self):
        result = list(wrap.cmalign(self.sequences))
        self.assertEqual(len(self.sequences), len(result))

    def test_mpi(self):
        result = list(wrap.cmalign(self.sequences, mpi_args=['-np', '2']))
        self.assertEqual(len(self.sequences), len(result))

class CMTestCase(unittest.TestCase):
    def test_find_cm(self):
        self.assertTrue(os.path.isfile(wrap.CM))

@unittest.skipUnless(util.which('FastTree'), "FastTree not found")
class AsRefpkgTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = SeqIO.parse(
                util.data_path('taxon171552.seqs.fasta'), 'fasta')

    def test_as_refpkg(self):
        with wrap.as_refpkg(self.sequences) as refpkg:
            self.assertTrue(os.path.isdir(refpkg.path))

            if util.which('rppr'):
                out = subprocess.check_output(['rppr', 'check', '-c', refpkg.path])
                self.assertTrue('OK!' in out, out)

def suite():
    s = unittest.TestSuite()
    classes = [CmAlignTestCase, CMTestCase]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))

    return s
