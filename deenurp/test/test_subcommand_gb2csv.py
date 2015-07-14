"""
Test gb2csv output, records.csv and references.csv
"""

import os
import unittest

from cStringIO import StringIO

from deenurp.subcommands import gb2csv


def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)


class Gb2csvTestCase(unittest.TestCase):

    def test01(self):
        """
        test records output
        """

        test_input = data_path('records.gb')
        test_reference = data_path('records.csv')

        class args(object):
            infile = test_input
            out = StringIO()
            references_out = None

        gb2csv.action(args)

        self.assertEqual(open(test_reference).read(), args.out.getvalue())

    def test02(self):
        """
        test references output
        """

        test_input = data_path('records.gb')
        test_reference = data_path('references.csv')

        class args(object):
            infile = test_input
            out = StringIO()
            references_out = StringIO()

        gb2csv.action(args)

        self.assertEqual(open(test_reference).read(),
                         args.references_out.getvalue())

    def test03(self):
        """
        test null case
        """

        test_reference = data_path('records.csv')

        class args(object):
            infile = StringIO('')
            out = StringIO()
            references_out = None

        gb2csv.action(args)

        self.assertEqual(next(open(test_reference)), args.out.getvalue())
