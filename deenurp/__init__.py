import os
import sys

try:
    with open(os.path.join(os.path.dirname(__file__), 'data', 'ver')) as v:
        __version__ = v.read().strip().lstrip('v')
except Exception, e:
    sys.stderr.write(e)
    __version__ = ''
