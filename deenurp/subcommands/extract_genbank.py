"""
Convert a genbank record to csv format

Output:

Genbank
[id, name, description, version, accession, gi, tax_id,
 comment, date, source, keywords, organism]

Seq_References
[version, title, authors, comment,
 consrtm, journal, medline_id, pubmed_id]
"""

import csv
import itertools
import operator
import sys

from deenurp import util, seq_info


def build_parser(parser):
    parser.add_argument('infile',
                        default=sys.stdin,
                        nargs='?',
                        type=util.file_opener(mode='r'),
                        help="""path to genbank file [default: stdin]""")
    parser.add_argument('--references-out',
                        type=util.file_opener(mode='w'),
                        help="""output references""")
    parser.add_argument('--database', help='Path to taxonomy database')
    parser.add_argument('--out',
                        default=sys.stdout,
                        type=util.file_opener(mode='w'),
                        help="""output path [default: stdout]""")


def action(args):
    if args.database:
        seq_info.set_taxonomy(args.database)

    records = seq_info.parse(args.infile, format='gb')
    records = util.Counter(records, prefix='Record ')

    out = csv.DictWriter(args.out, fieldnames=seq_info.info_fieldnames)
    out.writeheader()

    if args.references_out:
        fieldnames = seq_info.ref_fieldnames
        out_refs = csv.DictWriter(args.references_out, fieldnames=fieldnames)
        out_refs.writeheader()

    records = sorted(records, key=operator.attrgetter('version'))
    records = itertools.groupby(records, key=operator.attrgetter('version'))
    records = (list(recs) for ver, recs in records)

    for recs in records:
        if args.references_out:
            # write reference file before deduplicating records
            has_refs = [rec for rec in recs if rec.has_references()]
            refs = set([ref for rec in has_refs for ref in rec.references()])
            out_refs.writerows(ref.ref_dict() for ref in refs)
        out.writerows(rec.record_dict() for rec in set(recs))
