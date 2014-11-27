import unittest

modules = [
    'test_outliers',
    'test_search',
    'test_subcommand_hrefpkg_build',
    'test_subcommand_filter_outliers',
    'test_util',
    'test_wrap',
]

def suite():
    s = unittest.TestSuite()
    for m in modules:
        module = __import__(__name__ + '.' + m, fromlist=m)
        s.addTests(module.suite())

    return s
