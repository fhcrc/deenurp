import unittest

from deenurp import uclust
from deenurp.test import util


class ParseUclustAsDfTestCase(unittest.TestCase):
    def setUp(self):
        self.infile = util.data_path('fusobacterium_nucleatum_refs.uc')

    def test01(self):
        df = uclust.parse_uclust_as_df(self.infile)
        # target_label always has a value for types S and H
        self.assertFalse(any(df[df['type'] != 'C']['target_label'].isnull()))
