"""
Choose optimal reference sequences
"""

import argparse
import contextlib
import shlex
import sqlite3
from romperroom import RefsetInternalFasta

RefsetInternalFasta.install()

from Bio import SeqIO

from .. import select, search, wrap

def unique_sequences(it):
    seen = set()
    for sequence in it:
        if sequence.id in seen:
            continue
        seen.add(sequence.id)
        yield sequence

def build_parser(p):
    p.add_argument('search_db', help="""Output of `deenurp search-sequences`""")
    p.add_argument('--threads', help="""Number of threads [default:
            %(default)d]""", type=int, default=6)
    p.add_argument('--mpi-args', type=shlex.split, default=[])
    p.add_argument('--refs-per-cluster', type=int, default=5, help="""Maximum
            references per cluster [default: %(default)d]""")
    p.add_argument('--cluster-candidates', type=int, default=30,
            help="""Maximum candidate refs per cluster [default: %(default)d]""")
    p.add_argument('--min-mass-prop', help="""Minimum proportion of total mass
            in a cluster to require before searching for references [default:
            %(default)f]""", type=float, default=-1.0)
    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))

def action(args):
    with wrap.tempcopy(args.search_db) as search_path:
        s = sqlite3.connect(search_path)
        with contextlib.closing(s):
            search_db = search.open_database(s)
            sequences = select.choose_references(search_db,
                    args.refs_per_cluster, candidates=args.cluster_candidates,
                    threads=args.threads, min_cluster_prop=args.min_mass_prop, mpi_args=args.mpi_args)
            with args.output as fp:
                SeqIO.write(unique_sequences(sequences), fp, 'fasta')
