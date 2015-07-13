import unittest

modules = [
    'test_outliers',
    'test_search',
    'test_subcommand_hrefpkg_build',
    'test_subcommand_filter_outliers',
    'test_util',
    'test_wrap',
    'test_subcommand_gb2csv'
]


def get_test_suites(module):
    for name in dir(module):
        obj = getattr(module, name)
        try:
            if issubclass(obj, unittest.TestCase):
                yield unittest.makeSuite(obj)
        except TypeError:
            pass


def suite():
    s = unittest.TestSuite()
    for m in modules:
        module = __import__(__name__ + '.' + m, fromlist=m)
        for suite in get_test_suites(module):
            s.addTests(suite)

    return s
