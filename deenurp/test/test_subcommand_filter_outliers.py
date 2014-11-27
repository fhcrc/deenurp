import unittest

from deenurp.subcommands import filter_outliers
from deenurp.test import util

class FilterOutliersFunctions(unittest.TestCase):
    def setUp(self):
        pass

    def test01(self):
        filename = util.data_path('test_input.allpairs.blast6out')
        filter_outliers.parse_usearch_allpairs(filename)


def suite():
    s = unittest.TestSuite()
    classes = [FilterOutliersFunctions]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))
    return s
