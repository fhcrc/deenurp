"""
Convert a genbank record to csv format (seq_info.csv and seq_refs.csv)

Output:

seq_info.csv
[id, name, description, version, accession, gi, tax_id,
 comment, date, source, keywords, organism]

seq_refs.vsv
[version, title, authors, comment,
 consrtm, journal, medline_id, pubmed_id]
"""

import csv
import itertools
import operator
import sys

from Bio import SeqIO
from deenurp import util, seq_info


def build_parser(parser):
    # input
    parser.add_argument('infile',
                        default=sys.stdin,
                        nargs='?',
                        type=util.file_opener(mode='r'),
                        help="""path to genbank file [default: stdin]""")

    # outputs
    parser.add_argument('fasta_out',
                        type=util.file_opener(mode='w'),
                        help='output fasta')
    parser.add_argument('seqinfo_out',
                        metavar='CSV',
                        type=util.file_opener(mode='w'),
                        help='output seq_info')
    parser.add_argument('--references-out',
                        type=util.file_opener(mode='w'),
                        help="""output references""")

    # optionals
    parser.add_argument('--database', help="""Path to taxonomy database
                                              for validating tax_ids""")


def action(args):
    if args.database:
        seq_info.set_taxonomy(args.database)

    records = seq_info.parse(args.infile, format='gb')
    records = util.Counter(records, prefix='Record ')

    seqinfo_out = csv.DictWriter(args.seqinfo_out,
                                 fieldnames=seq_info.info_fieldnames)
    seqinfo_out.writeheader()

    if args.references_out:
        fieldnames = seq_info.ref_fieldnames
        out_refs = csv.DictWriter(args.references_out, fieldnames=fieldnames)
        out_refs.writeheader()

    # check out seqmagick.transform.sort to avoid writing records to memory
    records = sorted(records, key=operator.attrgetter('version'))
    records = itertools.groupby(records, key=operator.attrgetter('version'))
    records = (list(recs) for ver, recs in records)

    fa = []

    for recs in records:
        # write reference file before taking `set' of records
        if args.references_out:
            has_refs = [rec for rec in recs if rec.has_references()]
            refs = set([ref for rec in has_refs for ref in rec.references()])
            out_refs.writerows(ref.ref_dict() for ref in refs)
        recs = set(recs)
        seqinfo_out.writerows(rec.record_dict() for rec in recs)
        fa.extend(recs)

    SeqIO.write(fa, args.fasta_out, 'fasta')
