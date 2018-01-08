import sys
import unittest

from .. import test

if __name__ == '__main__':

    result = unittest.TestResult()
    suite = test.suite()
    outcome = suite.run(result)
    if outcome.wasSuccessful():
        print('ok: ' + str(outcome))
    else:
        print('--> {} failures:'.format(len(outcome.failures)))
        for testcase, tb in outcome.failures:
            msg = str(testcase)
            print('=' * len(msg))
            print(msg + '\n')
            print(tb.strip())
            print('=' * len(msg))
        sys.exit(1)
