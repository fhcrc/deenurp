import os.path
import unittest

from deenurp import expand

DN = os.path.dirname(__file__)

def data_path(*args):
    return os.path.join(DN, 'data', *args)

class NodeTestCase(unittest.TestCase):
    def setUp(self):
        with open(data_path('test_taxtable.csv')) as fp:
            self.root = expand.Node.of_taxtable(fp)

    def test_index(self):
        self.assertEqual(self.root, self.root.get_node('1'))
        self.assertEqual('1', self.root.tax_id)

    def test_iter(self):
        self.assertEqual(356, sum(1 for i in self.root))
