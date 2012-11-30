import functools
import os
import os.path

"""path to test file"""
data_path = functools.partial(os.path.join, os.path.dirname(__file__), 'data')
