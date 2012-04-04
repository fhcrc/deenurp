"""
Add all sequences in a sequence file whose tax_id intersects with a taxtable.
"""

import argparse
import csv

from .. import wrap

def build_parser(p):
    p.add_argument('fasta_file', help="""Sequence file to augment reference
            package with""", metavar='<fasta_file>')
    p.add_argument('tax_map', help="""CSV file containing [seqid,tax_id] rows
            for <fasta_file>""", type=argparse.FileType('r'))
    p.add_argument('--header', help="""Tax map has a header [default:
            %(default)s]""", default=False, action='store_true')
    p.add_argument('taxtable', help="""Output of `taxit taxtable`""", type=argparse.FileType('r'))
    p.add_argument('outfile', type=argparse.FileType('w'))

def action(args):
    # Load all tax_ids
    with args.taxtable as fp:
        r = csv.DictReader(fp)
        tax_ids = frozenset(i['tax_id'] for i in r)

    with args.tax_map:
        tax_map = wrap.load_tax_maps([args.tax_map], args.header)

    sequence_ids = (k for k, v in tax_map.items() if v in tax_ids)

    # Fetch
    with args.outfile:
        wrap.esl_sfetch(args.fasta_file, sequence_ids, args.outfile)
