import unittest

from taxtastic.taxtable import TaxNode
from deenurp.subcommands import hrefpkg_build

class FindNodesTestCase(unittest.TestCase):
    def setUp(self):
        self.taxonomy = TaxNode(rank='root', name='root', tax_id='1')
        self.taxonomy.ranks = ['root', 'class', 'genus', 'species']
        g1 = TaxNode(rank='genus', name='g1', tax_id='2')
        self.g1 = g1
        g1.sequence_ids = set(['s1', 's2'])
        self.taxonomy.add_child(g1)
        g1.add_child(TaxNode(rank='species', name='s1', tax_id='s1'))
        g1.add_child(TaxNode(rank='species', name='s2', tax_id='s2'))

        g2 = TaxNode(rank='genus', name='g2', tax_id='3')
        self.taxonomy.add_child(g2)
        s3 = TaxNode(rank='species', name='s3', tax_id='s3')
        s3.sequence_ids = set(['s3', 's4'])
        g2.add_child(s3)
        g2.add_child(TaxNode(rank='species', name='s4', tax_id='s4'))

    def test_find_nodes(self):
        r = list(hrefpkg_build.find_nodes(self.taxonomy, 'class'))
        self.assertEqual(frozenset(['g1', 's3']), frozenset(i.name for i in r))

    def test_find_nodes_below_rank(self):
        r = list(hrefpkg_build.find_nodes(self.taxonomy, 'genus'))
        self.assertEqual(['s3'], [i.name for i in r])

