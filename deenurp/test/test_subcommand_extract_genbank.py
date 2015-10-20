"""
Test extract-genbank output, records.csv and references.csv

FIXME:
"""

import os
import unittest

from cStringIO import StringIO

from deenurp.subcommands import ncbi_extract_genbank


def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)


class ExtractGenbankTestCase(unittest.TestCase):

    def test01(self):
        """
        test records output
        """

        test_input = data_path('records.gb')
        test_seqinfo = data_path('records.csv')
        test_fasta = data_path('records.fa')

        class args(object):
            infile = test_input
            seqinfo_out = StringIO()
            fasta_out = StringIO()
            references_out = None
            database = None

        ncbi_extract_genbank.action(args)

        self.assertEqual(open(test_fasta).read(),
                         args.fasta_out.getvalue())
        self.assertEqual(open(test_seqinfo).read(),
                         args.seqinfo_out.getvalue())

    def test02(self):
        """
        test references output
        """

        test_input = data_path('records.gb')
        test_reference = data_path('references.csv')

        class args(object):
            infile = test_input
            fasta_out = StringIO()
            seqinfo_out = StringIO()
            references_out = StringIO()
            database = None

        ncbi_extract_genbank.action(args)

        self.assertEqual(open(test_reference).read(),
                         args.references_out.getvalue())

    def test03(self):
        """
        test null case
        """

        test_seqinfo = data_path('records.csv')
        test_reference = data_path('references.csv')

        class args(object):
            infile = StringIO('')
            seqinfo_out = StringIO()
            fasta_out = StringIO()
            references_out = StringIO()
            database = None

        ncbi_extract_genbank.action(args)

        self.assertEqual(next(open(test_seqinfo)),
                         args.seqinfo_out.getvalue())
        self.assertEqual(next(open(test_reference)),
                         args.references_out.getvalue())
