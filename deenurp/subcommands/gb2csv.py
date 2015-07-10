"""
Convert a genbank record to csv format

Output:

Genbank
[id, name, description, version, accession, gi, tax_id,
 comment, date, source, keywords, organism]

References
[version, title, authors, comment,
 consrtm, journal, medline_id, pubmed_id]
"""

import argparse
import csv
import sys

from Bio import SeqIO

from deenurp.util import accession_version_of_genbank, tax_of_genbank


def build_parser(parser):
    parser.add_argument('infile',
                        default=sys.stdin,
                        nargs='?',
                        help="""path to genbank file [default: stdin]""")
    parser.add_argument('--references-out',
                        type=argparse.FileType('w'),
                        help="""output references""")
    parser.add_argument('--out',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="""output path [default: stdout]""")


def action(args):
    records = SeqIO.parse(args.infile, format='gb')

    fieldnames = ['version', 'accession', 'id', 'name', 'description',
                  'gi', 'tax_id', 'date', 'source', 'keywords', 'organism']

    out = csv.DictWriter(args.out, fieldnames=fieldnames)
    out.writeheader()

    if args.references:
        fieldnames = ['version', 'title', 'authors', 'comment',
                      'consrtm', 'journal', 'medline_id', 'pubmed_id']
        out_refs = csv.DictWriter(args.references_out, fieldnames=fieldnames)
        out_refs.writeheader()

    for record in records:
        accession, version = accession_version_of_genbank(record)
        annotations = record.annotations
        keywords = ';'.join(record.annotations.get('keywords', []))
        out.writerow(dict(accession=annotations['accessions'][0],
                          date=annotations['date'],
                          description=record.description,
                          gi=annotations.get('gi', ''),
                          id=record.id,
                          keywords=keywords,
                          name=record.name,
                          organism=annotations['organism'],
                          source=annotations['source'],
                          tax_id=tax_of_genbank(record),
                          version=version,
                          ))

        if out_refs:
            for ref in annotations['references']:
                out_refs.writerow(dict(authors=ref.authors,
                                       comment=ref.comment,
                                       consrtm=ref.consrtm,
                                       journal=ref.journal,
                                       medline_id=ref.medline_id,
                                       pubmed_id=ref.pubmed_id,
                                       title=ref.title,
                                       version=version))