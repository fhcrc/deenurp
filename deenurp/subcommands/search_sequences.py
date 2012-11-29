"""
Search a set of sequences against a sequence database for reference package
candidates.
"""

import argparse
import collections
import sqlite3

from .. import search

def build_parser(p):
    p.add_argument('sequence_file', help="""Fasta file containing query
            sequences""", metavar='<query_fasta>')
    p.add_argument('output', help="""Output database to write""", metavar='<output_db>')
    p.add_argument('ref_database', help="""Reference sequence database""")
    p.add_argument('ref_meta', help="""Reference sequence metadata""")
    p.add_argument('--weights', help="""Weights, in a `guppy dedup
            -m`-compatible dedup file""", type=argparse.FileType('r'))
    p.add_argument('--group-field', help="""Column to indicate group
            membership for a reference sequence (e.g., OTU; NCBI taxon id)
            [default: %(default)s]""", default='cluster')
    p.add_argument('--specimen-map',
            help="""CSV file containing two-column rows, with read name in the
            first, specimen identifier in the second [compatible with pplacer
            split placefiles]""", type=argparse.FileType('r'))
    uc = p.add_argument_group('UCLUST')
    uc.add_argument('--maxaccepts', default=5, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--maxrejects', default=40, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--search-identity', default=0.97, type=float,
            help="""Clustering identity level [default: %(default).2f]""")
    uc.add_argument('--select-threshold', help="""Select hits within
            %(metavar)s of best hit pct_id [default: %(default).2f]""",
            default=search.SELECT_THRESHOLD, metavar='THRESHOLD')

def action(args):
    con = sqlite3.connect(args.output)

    specimens = collections.defaultdict(str)
    if args.specimen_map:
        with args.specimen_map as fp:
            specimens = search.load_specimen_map(fp)
    else:
        specimens = collections.defaultdict(str)
    weights = None
    if args.weights:
        with args.weights:
            weights = search.dedup_info_to_counts(args.weights, specimens)
        assert weights
    search.create_database(con, args.sequence_file, ref_fasta=args.ref_database,
            ref_meta=args.ref_meta,
            weights=weights, maxaccepts=args.maxaccepts,
            maxrejects=args.maxrejects, search_id=args.search_identity,
            quiet=args.verbosity == 0, select_threshold=args.select_threshold,
            group_field=args.group_field)
