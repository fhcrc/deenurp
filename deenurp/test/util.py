import functools
import os
import os.path

from deenurp.util import which

"""path to test file"""
data_path = functools.partial(os.path.join, os.path.dirname(__file__), 'data')
