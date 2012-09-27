import unittest

modules = ['test_search', 'test_subcommand_hrefpkg_build', 'test_util', 'test_wrap']

def suite():
    s = unittest.TestSuite()
    for m in modules:
        module = __import__(__name__ + '.' + m, fromlist=m)
        s.addTests(module.suite())

    return s
