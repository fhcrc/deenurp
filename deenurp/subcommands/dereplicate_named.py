"""Dereplicate named reference sequences prior to filtering

"""

import csv
import logging
import shutil

import pandas as pd

from Bio import SeqIO
from deenurp import uclust
from taxtastic.taxtable import TaxNode

from .. import wrap, util


def build_parser(parser):
    parser.add_argument('seqs', help="""Named sequences""")
    parser.add_argument('seq_info', help="""Sequence info file""")
    parser.add_argument('-t', '--taxonomy',
                        help="""Taxonomy as taxtable; optional
                        if a grouping term is available in seq_info""")
    parser.add_argument('-g', '--group-on', default='species',
                        help="""Field in seq_info on which to group sequences""")
    parser.add_argument('-i', '--cluster-id', default=0.985,
                        type=float, help="""Cluster ID [default: %(default).3f]""")
    parser.add_argument('--seq-info-out')


def action(args):

    dtype = {'gi': str, 'tax_id': str, 'species': str}
    seq_info = pd.read_csv(args.seq_info, dtype=dtype, index_col='seqname')

    # TODO: join with taxonomy if provided
    grouped = seq_info.groupby(args.group_on, sort=False)

    for key, grp in grouped:
        print key, grp.shape
