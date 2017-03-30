"""
Choose reference sequences.
"""

import argparse
import contextlib
import csv
import operator
import sqlite3

from Bio import SeqIO

from .. import config, search, select, util


def meta_writer(fp):
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(('seqname', 'cluster_id', 'max_weight', 'mean_weight'))

    def inner(sequences):
        with fp:
            for sequence in sequences:
                writer.writerow(
                    (sequence.id,
                     sequence.annotations['cluster_name'],
                     sequence.annotations['max_weight'],
                     sequence.annotations['mean_weight']))
                yield sequence

    return inner


class track_attr(object):
    """
    Track every unique value seen in ``attr`` in an iterable

    Results are stored in ``seen``
    """

    def __init__(self, attr, iterable):
        self.attr = attr
        self.iterable = iterable
        self.seen = set()

    def __iter__(self):
        for i in self.iterable:
            self.seen.add(getattr(i, self.attr))
            yield i


def build_parser(p):
    p.add_argument(
        'search_db',
        help="""Output of `deenurp search-sequences`""")
    p.add_argument(
        'output',
        type=argparse.FileType('w'),
        help="Output file (fasta)")

    p.add_argument(
        '--threads',
        help="""Number of threads [default:%(default)d]""",
        type=int,
        default=config.DEFAULT_THREADS)

    selection_options = p.add_argument_group('Selection Options')
    selection_options.add_argument(
        '--refs-per-cluster', metavar='INT',
        type=int,
        default=5,
        help="""Maximum references per cluster [default: %(default)d]""")
    selection_options.add_argument(
        '--min-mass-prop', metavar='FLOAT',
        help="""Minimum proportion of total mass in
        a cluster to require before including references [default:
        %(default)f]""",
        type=float,
        default=select.MIN_CLUSTER_PROP)
    selection_options.add_argument(
        '--include-clusters', metavar='FILE',
        type=argparse.FileType('r'),
        help="""Select sequences for cluster IDs in %(metavar)s,
        regardless of whether they had hits among the query
        sequences""")
    selection_options.add_argument(
        '--exclude-clusters', metavar='FILE',
        type=argparse.FileType('r'),
        help=('List of cluster identifiers to exclude from the results'))
    # selection_options.add_argument(
    #     '--include-sequences',
    #     type=argparse.FileType('r'),
    #     help=('List of sequences to include in the results [NOT IMPLEMENTED]'))
    selection_options.add_argument(
        '--exclude-sequences', metavar='FILE',
        type=argparse.FileType('r'),
        help=('List of sequence ids to exclude from the results'))

    info_options = p.add_argument_group('Sequence info options')
    info_options.add_argument(
        '--seqinfo-out', metavar='FILE',
        type=argparse.FileType('w'),
        help="""File to write merged metadata""")
    info_options.add_argument(
        '--output-meta', metavar='FILE',
        help="""File to write selection metadata""",
        type=argparse.FileType('w'))


def extract_meta(ids, search_db, out_fp):
    """
    Subset sequence metadata to ids, writing the results to out_fp
    """
    params = search.load_params(search_db)
    seqinfo = params['ref_meta']
    with open(seqinfo) as fp:
        r = csv.DictReader(fp)
        rows = (i for i in r if i['seqname'] in ids)
        w = csv.DictWriter(out_fp, r.fieldnames, lineterminator='\n',
                           quoting=csv.QUOTE_NONNUMERIC)
        w.writeheader()
        w.writerows(rows)


def action(args):
    include_clusters = None
    if args.include_clusters:
        with args.include_clusters as fp:
            include_clusters = set(fp.read().split())

    exclude_clusters = None
    if args.exclude_clusters:
        with args.exclude_clusters as fp:
            exclude_clusters = set(fp.read().split())

    # include_sequences = None
    # if args.include_sequences:
    #     with args.include_sequences as fp:
    #         include_sequences = set(fp.read().split())

    exclude_sequences = None
    if args.exclude_sequences:
        with args.exclude_sequences as fp:
            exclude_sequences = set(fp.read().split())

    with util.tempcopy(args.search_db) as search_path:
        search_db = sqlite3.connect(search_path)
        with contextlib.closing(search_db):
            sequences = select.choose_references(
                search_db,
                args.refs_per_cluster,
                threads=args.threads,
                min_cluster_prop=args.min_mass_prop,
                include_clusters=include_clusters,
                exclude_clusters=exclude_clusters,
                # include_sequences=include_sequences,
                exclude_sequences=exclude_sequences)

            with args.output as fp:
                # Unique IDs
                sequences = util.unique(
                    sequences, key=operator.attrgetter('id'))
                sequences = util.unique(sequences, key=lambda s: str(s.seq))
                if args.output_meta:
                    sequences = meta_writer(args.output_meta)(sequences)
                sequences = track_attr('id', sequences)
                SeqIO.write(sequences, fp, 'fasta')

            if args.seqinfo_out:
                seen_ids = sequences.seen
                extract_meta(seen_ids, search_db, args.seqinfo_out)
