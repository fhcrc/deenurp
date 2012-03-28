import argparse
import contextlib
import sqlite3

from Bio import SeqIO

from .. import select, search

def build_parser(p):
    p.add_argument('search_db', help="""Output of `deenurp search-sequences`""",
            type=sqlite3.connect)
    p.add_argument('ref_db', help="""Reference sequence database""", type=sqlite3.connect)
    p.add_argument('--threads', help="""Number of threads [default:
            %(default)d]""", type=int, default=6)
    p.add_argument('--refs-per-cluster', type=int, default=5, help="""Maximum
            references per cluster [default: %(default)d]""")
    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))

def action(args):
    search_db = search.open_database(args.search_db)
    with contextlib.closing(args.search_db), contextlib.closing(args.ref_db):
        sequences = select.choose_references(search_db, args.ref_db,
                args.refs_per_cluster, args.threads)
        with args.output as fp:
            SeqIO.write(sequences, fp, 'fasta')
