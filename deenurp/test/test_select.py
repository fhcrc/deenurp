import os.path
import unittest

from Bio import SeqIO

from deenurp import select

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

