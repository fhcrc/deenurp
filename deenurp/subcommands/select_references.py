"""
Choose reference sequences
"""

import argparse
import contextlib
import csv
import shlex
import sqlite3
from romperroom import RefsetInternalFasta

RefsetInternalFasta.install()

from Bio import SeqIO

from .. import select, search, wrap

def _load_tax_maps(fps, has_header=False):
    d = {}
    for fp in fps:
        reader = csv.reader(fp)
        if has_header:
            next(reader) # Skip
        for row in reader:
            name, taxid = row[:2]
            if name in d and taxid != d[name]:
                raise ValueError("Multiple tax_ids specified for {0}".format(name))
            d[name] = taxid
    return d

def unique_sequences(it):
    seen = set()
    for sequence in it:
        if sequence.id in seen:
            continue
        seen.add(sequence.id)
        yield sequence

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

    taxonomy_options = p.add_argument_group('Taxonomy Options')
    taxonomy_options.add_argument('-m', '--taxid-map', help="""CSV File mapping
            from sequence ID to tax_id. May be specified multiple times.""",
            action='append', dest='taxid_maps', metavar='MAP', type=argparse.FileType('r'))
    taxonomy_options.add_argument('--map-header', help="""Taxonomy map(s) have
            a header row. [default: %(default)s]""", default=False,
            action='store_true')
    taxonomy_options.add_argument('--tax-map-out', type=argparse.FileType('w'),
            help="""File to write output seqid,taxid map""")

def action(args):
    taxid_map = None
    if args.taxid_maps and not args.tax_map_out:
        raise argparse.ArgumentError("--tax-map-out required with --taxid-map")
    if args.taxid_maps:
        taxid_map = _load_tax_maps(args.taxid_maps, args.map_header)
        headers = ('seqname', 'tax_id')
        writer = csv.writer(args.tax_map_out, lineterminator='\n')
        writer.writerow(headers)

        def write_taxid(sequences):
            with args.tax_map_out:
                for sequence in sequences:
                    writer.writerow((sequence.id, taxid_map.get(sequence.id)))
                    yield sequence
    with wrap.tempcopy(args.search_db) as search_path:
        s = sqlite3.connect(search_path)
        with contextlib.closing(s):
            search_db = search.open_database(s)
            sequences = select.choose_references(search_db,
                    args.refs_per_cluster, candidates=args.cluster_candidates,
                    threads=args.threads, min_cluster_prop=args.min_mass_prop,
                    mpi_args=args.mpi_args, cluster_factor=args.cluster_factor)
            with args.output as fp:
                sequences = unique_sequences(sequences)
                if taxid_map:
                    sequences= write_taxid(sequences)
                SeqIO.write(sequences, fp, 'fasta')
