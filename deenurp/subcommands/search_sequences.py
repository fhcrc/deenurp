"""
Search a set of sequences against a sequence database for reference package
candidates.
"""

import argparse
import sqlite3

from .. import search

def build_parser(p):
    p.add_argument('sequence_file', help="""Fasta file containing query
            sequences""", metavar='<query_fasta>')
    p.add_argument('output', help="""Output database to write""", metavar='<output_db>')
    p.add_argument('sequence_databases', metavar='<db_fasta>', help="""Sequence
            database for possible reference sequences. Databases are searched
            in order. Any clusters with hits to an earlier database are not
            re-searched.""", nargs='+')
    p.add_argument('--weights', help="""Weights, in a `guppy dedup
            -m`-compatible dedup file""", type=argparse.FileType('r'))
    uc = p.add_argument_group('UCLUST')
    uc.add_argument('--maxaccepts', default=5, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--maxrejects', default=40, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--cluster-identity', default=0.99, type=float,
            help="""Clustering identity level [default: %(default).2f]""")
    uc.add_argument('--search-identity', default=0.99, type=float,
            help="""Clustering identity level [default: %(default).2f]""")

def action(args):
    con = sqlite3.connect(args.output)
    weights = None
    if args.weights:
        with args.weights:
            weights = search.dedup_info_to_counts(args.weights)
        assert weights
    search.create_database(con, args.sequence_file, args.sequence_databases,
            weights=weights, maxaccepts=args.maxaccepts,
            maxrejects=args.maxrejects, cluster_id=args.cluster_identity,
            search_id=args.search_identity, quiet=args.verbosity == 0)
