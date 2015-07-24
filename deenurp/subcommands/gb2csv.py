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


class Seq_Info(SeqRecord.SeqRecord):
    attrs = ['version', 'accession', 'id', 'name', 'description',
             'gi', 'tax_id', 'date', 'source', 'keywords', 'organism']

    def rec_dict(self):
        return dict((a, self.__getattribute__(a)) for a in Seq_Info.attrs)

    def __eq__(self, other):
        return isinstance(other, Seq_Ref) and self.hash() == other.hash()

    def __hash__(self):
        """
        Unique Seq_Info ID

        Uniqueness limited to id, version and sequence
        """
        if not hasattr(self, 'hash_id'):
            self.hash_id = hash((self.id, self.version))
        return self.hash_id

    def hash(self):
        return self.__hash__()

    def has_references(self):
        return 'references' in self.annotations

    def references(self):
        """
        Might want to cache this in a variable if referencing more than once
        """
        if self.has_references():
            refs = self.annotations['references']
            return [convert_record(ref, Seq_Ref, self.version) for ref in refs]
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


class Seq_Ref(SeqFeature.Reference):
    attrs = ['version', 'title', 'authors', 'comment',
             'consrtm', 'journal', 'medline_id', 'pubmed_id']

    def __hash__(self):
        """
        Unique Seq_Ref ID

        Combo of all ref_names attributes
        """
        if not hasattr(self, 'hash_id'):
            attributes = tuple(self.__getattribute__(a) for a in Seq_Ref.attrs)
            self.hash_id = hash(attributes)
        return self.hash_id

    def __eq__(self, other):
        return isinstance(other, Seq_Ref) and self.hash() == other.hash()

    def hash(self):
        return self.__hash__()

    def ref_dict(self):
        return dict((a, self.__getattribute__(a)) for a in Seq_Ref.attrs)

    def setAttributes(self, version):
        title = self.title
        self.title = title if re.search('\w', title) else ''
        self.version = version
        return self


def convert_record(record, cls, *attrs):
    record.__class__ = cls
    return record.setAttributes(*attrs)


def action(args):
    records = SeqIO.parse(args.infile, format='gb')
    records = (convert_record(record, Seq_Info) for record in records)
    records = Counter(records, prefix='Record ')

    out = csv.DictWriter(args.out, fieldnames=Seq_Info.attrs)
    out.writeheader()

    if args.references_out:
        fieldnames = Seq_Ref.attrs
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
        out.writerows(rec.rec_dict() for rec in set(recs))
