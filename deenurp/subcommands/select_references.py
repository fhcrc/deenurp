"""
Choose optimal reference sequences
"""

import argparse
import contextlib
import sqlite3
from romperroom import RefsetInternalFasta

RefsetInternalFasta.install()

from Bio import SeqIO

from .. import select, search, wrap

def build_parser(p):
    p.add_argument('search_db', help="""Output of `deenurp search-sequences`""")
    p.add_argument('ref_db', help="""Reference sequence database""")
    p.add_argument('--threads', help="""Number of threads [default:
            %(default)d]""", type=int, default=6)
    p.add_argument('--refs-per-cluster', type=int, default=5, help="""Maximum
            references per cluster [default: %(default)d]""")
    p.add_argument('--cluster-candidates', type=int, default=30,
            help="""Maximum candidate refs per cluster [default:
            %(default)d]""")
    p.add_argument('--min-mass-prop', help="""Minimum proportion of total mass
            in a cluster to require before searching for references [default:
            %(default)f]""", type=float, default=0.0)
    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))

def action(args):
    with wrap.tempcopy(args.search_db) as search_path, wrap.tempcopy(args.ref_db) as ref_path:
        s = sqlite3.connect(search_path)
        ref = sqlite3.connect(ref_path)
        with contextlib.closing(s), contextlib.closing(ref):
            search_db = search.open_database(s)
            sequences = select.choose_references(search_db, ref,
                    args.refs_per_cluster, candidates=args.cluster_candidates,
                    threads=args.threads, min_cluster_prop=args.min_mass_prop)
            with args.output as fp:
                SeqIO.write(sequences, fp, 'deenurp-internal-fasta')
