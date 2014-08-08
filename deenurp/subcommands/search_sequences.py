"""Search a set of sequences against a sequence database for reference package
candidates.

The search algorithm is as follows:

* A database search is performed against ``ref_database`` using
  ``uclust --pct_id=SEARCH_THRESHOLD``.
* Hits with pairwise identity less than ``SEARCH_IDENTITY`` are also
  discarded (note that search sensitivity can be improved by choosing
  a SEARCH_THRESHOLD < SEARCH_IDENTITY).
* For each query sequence, the hit with the highest pairwise identity
  is identified. The difference between the highest identity and the
  pairwise identity score of each additional hit is calculated, and
  those with a difference > ``SELECT_THRESHOLD`` are also discarded.
* For each query sequence, each group identified by ``--group-field``
  is represented by a single reference sequence only.

"""

import argparse
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
    p.add_argument('--sample-map',
            help="""CSV file containing two-column rows, with read name in the
            first, sample identifier in the second [compatible with pplacer
            split placefiles]""", type=argparse.FileType('r'))
    p.add_argument('--blacklist', type=argparse.FileType('r'),
            help="""List of cluster identifiers not to include in the results""")
    uc = p.add_argument_group('UCLUST')
    uc.add_argument('--maxaccepts', default=5, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--maxrejects', default=40, type=int,
            help="""[default: %(default)d]""")
    uc.add_argument('--search-threshold', help="""
            Minimum threshold for database search (ie, "uclust --id") [default: %(default).2f]""",
            default=search.SEARCH_THRESHOLD, metavar='THRESHOLD')
    uc.add_argument('--search-identity', default=search.SEARCH_IDENTITY, type=float,
            help="""Identity threshold for filtering results of database search
            [default: %(default).2f]""")
    uc.add_argument('--select-threshold', help="""Select hits within
            %(metavar)s of best hit pct_id [default: %(default).2f]""",
            default=search.SELECT_THRESHOLD, metavar='THRESHOLD')


def action(args):
    con = sqlite3.connect(args.output)

    blacklist = set()
    if args.blacklist:
        with args.blacklist:
            blacklist = set(i.strip() for i in args.blacklist)

    samples = None
    if args.sample_map:
        with args.sample_map as fp:
            samples = search.load_sample_map(fp)
    weights = None
    if args.weights:
        with args.weights:
            weights = search.dedup_info_to_counts(args.weights, samples)
        assert weights

    # create_database(search_threshold=args.search_threshold)
    # --> _search(search_threshold=search_threshold)
    # --> uclust.search(pct_id=search_threshold)
    # --> "uclust --id pct_id" and filter records after search

    # create_database(search_id=args.search_identity)
    # --> _create_tables(search_id=search_id)
    # --> (saved in params)
    # --> used to filter output of uclust_search in _search()

    # create_database(select_threshold=args.select_threshold)
    # --> _search(select_threshold=select_threshold) -->
    # select_hits(threshold=select_threshold)

    search.create_database(
        con,
        args.sequence_file,
        ref_fasta=args.ref_database,
        ref_meta=args.ref_meta,
        weights=weights,
        maxaccepts=args.maxaccepts,
        maxrejects=args.maxrejects,
        search_id=args.search_identity,
        quiet=args.verbosity == 0,
        select_threshold=args.select_threshold,
        search_threshold=args.search_threshold,
        group_field=args.group_field,
        blacklist=blacklist)
