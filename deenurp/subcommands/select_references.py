"""
Choose reference sequences.
"""

import argparse
import contextlib
import csv
import operator
import sqlite3

from Bio import SeqIO

from .. import search, select, util

def meta_writer(fp):
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(('seqname', 'cluster_id', 'max_weight', 'mean_weight'))

    def inner(sequences):
        with fp:
            for sequence in sequences:
                writer.writerow((sequence.id,
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
    p.add_argument('search_db', help="""Output of `deenurp search-sequences`""")
    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))

    p.add_argument('--threads', help="""Number of threads [default:
            %(default)d]""", type=int, default=6)

    selection_options = p.add_argument_group('Selection Options')
    selection_options.add_argument('--refs-per-cluster', type=int, default=5,
            help="""Maximum references per cluster [default: %(default)d]""")
    selection_options.add_argument('--min-mass-prop', help="""Minimum
            proportion of total mass in a cluster to require before including
            references [default: %(default)f]""", type=float, default=-1.0)

    info_options = p.add_argument_group('Sequence info options')
    info_options.add_argument('--seqinfo-out', type=argparse.FileType('w'),
            help="""File to write merged metadata""")
    info_options.add_argument('--output-meta', help="""File to write selection metadata""",
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
    with util.tempcopy(args.search_db) as search_path:
        search_db = sqlite3.connect(search_path)
        with contextlib.closing(search_db):
            sequences = select.choose_references(search_db,
                    args.refs_per_cluster,
                    threads=args.threads, min_cluster_prop=args.min_mass_prop)

            with args.output as fp:
                # Unique IDs
                sequences = util.unique(sequences, key=operator.attrgetter('id'))
                sequences = util.unique(sequences, key=lambda s: str(s.seq))
                if args.output_meta:
                    sequences = meta_writer(args.output_meta)(sequences)
                sequences = track_attr('id', sequences)
                SeqIO.write(sequences, fp, 'fasta')

            if args.seqinfo_out:
                seen_ids = sequences.seen
                extract_meta(seen_ids, search_db, args.seqinfo_out)
