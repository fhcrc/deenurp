import os.path
import operator
import unittest

from deenurp import util

class UniqueTestCase(unittest.TestCase):
    def test_nokey(self):
        l = [1, 2, 3, 1, 2, 4, 6, 7, 5, 1]
        expected = [1, 2, 3, 4, 6, 7, 5]
        actual = list(util.unique(l))
        self.assertEqual(expected, actual)

    def test_key(self):
        keys = ('n', 's')
        v = [(1, 'test'), (2, 'test'), (2, 'other')]
        l = [dict(zip(keys, i)) for i in v]
        expected1 = [{'n': 1, 's': 'test'},
                    {'n': 2, 's': 'test'}]
        actual1 = util.unique(l, key=operator.itemgetter('n'))
        self.assertEqual(expected1, list(actual1))

class MemoizeTestCase(unittest.TestCase):
    def test_function(self):
        d = {'test': object()}
        m = util.memoize(d.get)

        self.assertEqual(d['test'], m('test'))

        # Check memoize
        expected = d.pop('test')

        self.assertEqual(expected, m('test'))

class MaybeTempFileTestCase(unittest.TestCase):
    def test_tempfile(self):
        with util.maybe_tempfile(prefix='tmp') as tf:
            self.assertFalse(tf.closed)
            self.assertTrue(os.path.basename(tf.name).startswith('tmp'))
        self.assertTrue(tf.closed)
        self.assertFalse(os.path.exists(tf.name))

    def test_obj(self):
        o = object()
        with util.maybe_tempfile(o, prefix='tmp') as tf:
            self.assertEqual(o, tf)

def suite():
    s = unittest.TestSuite()
    classes = [MaybeTempFileTestCase, MemoizeTestCase, UniqueTestCase]
    for cls in classes:
        s.addTests(unittest.makeSuite(cls))

    return s
