"""
Choose reference sequences.
"""

import argparse
import contextlib
import csv
import operator
import shlex
import sqlite3

from Bio import SeqIO

from .. import select, search, wrap

def meta_writer(fp):
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(('seqname', 'cluster_id', 'weight_prop'))

    def inner(sequences):
        with fp:
            for sequence in sequences:
                writer.writerow((sequence.id,
                    sequence.annotations['cluster_id'],
                    sequence.annotations['weight_prop']))
                yield sequence

    return inner

def build_parser(p):
    p.add_argument('search_db', help="""Output of `deenurp search-sequences`""")
    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))

    mpi_group = p.add_argument_group('Number of processors')
    mpi_group.add_argument('--threads', help="""Number of threads [default:
            %(default)d]""", type=int, default=6)
    mpi_group.add_argument('--mpi-args', type=shlex.split, default=[])

    selection_options = p.add_argument_group('Selection Options')
    selection_options.add_argument('--refs-per-cluster', type=int, default=5,
            help="""Maximum references per cluster [default: %(default)d]""")
    selection_options.add_argument('--cluster-candidates', type=int,
            default=30, help="""Maximum candidate refs per cluster [default:
            %(default)d]""")
    selection_options.add_argument('--cluster-factor', type=int, default=2,
            help="""Factor by which to expand reference selection for merged
            clusters [default: %(default)d]""")
    selection_options.add_argument('--min-mass-prop', help="""Minimum
            proportion of total mass in a cluster to require before searching
            for references [default: %(default)f]""", type=float, default=-1.0)

    info_options = p.add_argument_group('Sequence info options')
    info_options.add_argument('--seqinfo-out', type=argparse.FileType('w'),
            help="""File to write merged metadata""")
    info_options.add_argument('--output-meta', help="""File to write selection metadata""",
            type=argparse.FileType('w'))

def action(args):
    with wrap.tempcopy(args.search_db) as search_path:
        s = sqlite3.connect(search_path)
        with contextlib.closing(s):
            search_db = search.open_database(s)
            sequences = select.choose_references(search_db,
                    args.refs_per_cluster, candidates=args.cluster_candidates,
                    threads=args.threads, min_cluster_prop=args.min_mass_prop,
                    mpi_args=args.mpi_args, cluster_factor=args.cluster_factor)

            with args.output as fp:
                # Unique IDs
                sequences = wrap.unique(sequences, key=operator.attrgetter('id'))
                sequences = wrap.unique(sequences, key=lambda s: str(s.seq))
                if args.output_meta:
                    sequences = meta_writer(args.output_meta)(sequences)
                SeqIO.write(sequences, fp, 'fasta')

            if args.seqinfo_out:
                select.merge_meta(fp.name, search_db, args.seqinfo_out)
