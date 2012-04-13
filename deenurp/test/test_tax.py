import os.path
import unittest

from deenurp import tax2tree

DN = os.path.dirname(__file__)

def data_path(*args):
    return os.path.join(DN, 'data', *args)

class NodeTestCase(unittest.TestCase):
    def setUp(self):
        with open(data_path('test_taxtable.csv')) as fp:
            self.root = tax2tree.TaxNode.from_taxtable(fp)

    def test_index(self):
        self.assertEqual(self.root, self.root.get_node('1'))
        self.assertEqual('1', self.root.tax_id)

    def test_iter(self):
        self.assertEqual(356, sum(1 for i in self.root))

    def test_lineage(self):
        node = self.root.get_node('1303')
        lineage = node.lineage()
        self.assertEqual(['1', '131567', '2', '1239', '91061', '186826', '1300', '1301', '1303'],
                [i.tax_id for i in lineage])
