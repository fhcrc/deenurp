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

import csv
import itertools
import re
import sys

from Bio import SeqIO, SeqRecord, SeqFeature

from deenurp.util import (accession_version_of_genbank,
                          tax_of_genbank, Counter, file_opener)


def build_parser(parser):
    parser.add_argument('infile',
                        default=sys.stdin,
                        nargs='?',
                        type=file_opener(mode='r'),
                        help="""path to genbank file [default: stdin]""")
    parser.add_argument('--references-out',
                        type=file_opener(mode='w'),
                        help="""output references""")
    parser.add_argument('--out',
                        default=sys.stdout,
                        type=file_opener(mode='w'),
                        help="""output path [default: stdout]""")


class SeqRec(SeqRecord.SeqRecord):
    rec_names = ['version', 'accession', 'id', 'name', 'description',
                 'gi', 'tax_id', 'date', 'source', 'keywords', 'organism']

    def rec_dict(self):
        return dict((a, self.__getattribute__(a)) for a in SeqRec.rec_names)

    def __hash__(self):
        """
        Unique SeqRec ID

        Combo of id, version and sequence should define uniqueness
        """
        return hash((self.id, self.version, str(self.seq)))

    def hash(self):
        if not hasattr(self, 'hash_id'):
            self.hash_id = self.__hash__()
        return self.hash_id

    def references(self):
        if 'references' in self.annotations:
            refs = self.annotations['references']
            return set(convert_record(ref, Ref) for ref in refs)
        else:
            return None

    def setAttributes(self):
        self.accession, self.version = accession_version_of_genbank(self)
        self.date = self.annotations['date']
        self.gi = self.annotations.get('gi', '')
        self.keywords = ';'.join(self.annotations.get('keywords', []))
        self.organism = self.annotations['organism']
        self.source = self.annotations['source']
        self.tax_id = tax_of_genbank(self)
        return self


class Ref(SeqFeature.Reference):
    ref_names = ['version', 'title', 'authors', 'comment',
                 'consrtm', 'journal', 'medline_id', 'pubmed_id']

    def __hash__(self):
        """
        Unique Ref ID

        Combo of all ref_names attributes
        """
        return hash([self.__getattribute__(a) for a in SeqRec.ref_names])

    def hash(self):
        if not hasattr(self, 'hash_id'):
            self.hash_id = self.__hash__()
        return self.hash_id

    def ref_dict(self):
        return dict((a, self.__getattribute__(a)) for a in SeqRec.ref_names)

    def setAttributes(self, record_id):
        title = self.title
        self.title = title if re.search('\w', title) else ''
        self.id = record_id
        return self


def convert_record(record, cls):
    record.__class__ = cls
    return record.setAttributes()


def action(args):
    records = SeqIO.parse(args.infile, format='gb')
    records = (convert_record(record, SeqRec) for record in records)
    records = set(records)
    records = Counter(records, prefix='Record ')

    out = csv.DictWriter(args.out, fieldnames=SeqRec.re_names)
    out.writeheader()

    if args.references_out:
        fieldnames = SeqRec.rec_names
        out_refs = csv.DictWriter(args.references_out, fieldnames=fieldnames)
        out_refs.writeheader()

    records = sorted(records, key=records.version)
    records = itertools.groupby(records, by=records.version)

    for _, records in records:
        out.writerows(r.rec_dict for r in records)

        if args.references_out:
            out_refs.writerows(ref.ref_dict for rec in records for ref in rec)
