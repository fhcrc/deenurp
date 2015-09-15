"""
Seq_Info object for output of Genbank records and reference information
"""

import Bio
import logging
import re
import sqlalchemy

from sqlalchemy.sql import select
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

from deenurp import util

_taxonomy = None

info_fieldnames = ['seqname', 'version', 'accession', 'name', 'description',
                   'gi', 'tax_id', 'date', 'source', 'keywords', 'organism',
                   'length', 'ambig_count', 'is_type', 'taxid_classified']

ref_fieldnames = ['version', 'title', 'authors', 'comment',
                  'consrtm', 'journal', 'medline_id', 'pubmed_id']


def set_taxonomy(database):
    global _taxonomy
    engine = sqlalchemy.create_engine('sqlite:///{0}'.format(database))
    _taxonomy = Taxonomy(engine, ncbi.ranks)
    return


def is_type(record):
    """
    Returns a boolean indicating whether a sequence is a member of a type
    strain, as indicated by the presence of the string '(T)' within the
    record description.
    """
    type_keywords = ['(T)', 'ATCC', 'NCTC', 'NBRC', 'CCUG',
                     'DSM', 'JCM', 'NCDO', 'NCIB', 'CIP']

    for t in type_keywords:
        if t in record.description:
            return True

    return False


def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record, returning None if no taxonomy is
    available.
    """
    # Check for bad name
    try:
        source = next(i for i in gb.features if i.type == 'source')
        taxon = next(i[6:] for i in source.qualifiers.get('db_xref', [])
                     if i.startswith('taxon:'))
        return taxon
    except StopIteration:
        return


def accession_version_of_genbank(record):
    """
    Return the accession and version of a Bio.SeqRecord.SeqRecord
    """
    annotations = record.annotations
    accession = annotations.get('accessions', [''])[0]
    if accession:
        if 'sequence_version' in annotations:
            version = int(annotations.get('sequence_version'))
        elif record.id.startswith(accession + '.'):
            version = int(record.id.split('.', 1)[1])
        else:
            version = 1
        version = '{}.{}'.format(accession, version)
    else:
        version = ''
    return accession, version


@util.memoize
def fetch_tax_info(tax_id):
    if _taxonomy is not None:
        c = _taxonomy.nodes.c  # columns
        s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
        s = s.where(c.tax_id == tax_id)
        return s.execute().fetchone()
    else:
        return tax_id, True, None, None


@util.memoize
def species_is_classified(tax_id):
    """
    return is_valid from ncbi taxonomy
    """
    res = fetch_tax_info(tax_id)
    if not res:
        return False

    tax_id, is_valid, parent_id, rank = res
    if rank == 'species' or rank is None:
        return is_valid
    elif tax_id == parent_id:
        return False
    else:
        # instead of recursion create a table of ranks by index. code
        # will be something liek ranks.index(rank) > species_index
        return species_is_classified(parent_id)


def parse(handle, *args, **kwargs):
    records = Bio.SeqIO.parse(handle, *args, **kwargs)
    return (convert_record(r, Seq_Info) for r in records)


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def convert_record(record, cls, *attrs):
    record.__class__ = cls
    return record.setAttributes(*attrs)


@util.memoize
def update_taxid(tax_id, name):
    """
    Update the tax_id if it has been merged.  If the tax_id is not in the
    database a tax_id will be assigned based on the organism name
    """
    try:
        if _taxonomy is not None:
            _taxonomy._node(tax_id)
    except KeyError as err:
        new_tax_id = _taxonomy._get_merged(tax_id)
        if new_tax_id != tax_id:
            msg = 'tax_id {} merged to {}'.format(tax_id, new_tax_id)
            logging.warn(msg)
            tax_id = new_tax_id
        elif name:
            try:
                tax_id, _, _ = _taxonomy.primary_from_name(name)
            except KeyError as err:
                logging.warn(err)
                tax_id = None
        else:
            msg = 'taxid {} not found in taxonomy, dropping'
            logging.warn(msg.format(tax_id))
            tax_id = None

    return tax_id


class Seq_Info(Bio.SeqRecord.SeqRecord):
    def record_dict(self):
        return dict((a, self.__getattribute__(a)) for a in info_fieldnames)

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
        self.length = len(self)
        self.ambig_count = count_ambiguous(str(self.seq))
        self.is_type = is_type(self)
        self.date = self.annotations['date']
        self.gi = self.annotations.get('gi', '')
        self.keywords = ';'.join(self.annotations.get('keywords', []))
        self.organism = self.annotations['organism']
        self.source = self.annotations['source']
        self.tax_id = update_taxid(tax_of_genbank(self), self.organism)
        self.taxid_classified = species_is_classified(self.tax_id)
        return self

    def fasta_str(self):
        if not hasattr(self, self.fasta):
            self.fasta = '{name}\n{seq.seq}'.format(**self.__dict__)
        return self.fasta


class Seq_Ref(Bio.SeqFeature.Reference):
    def __hash__(self):
        """
        Unique Seq_Ref ID

        Combo of all ref_names attributes
        """
        if not hasattr(self, 'hash_id'):
            attributes = tuple(self.__getattribute__(a)
                               for a in ref_fieldnames)
            self.hash_id = hash(attributes)
        return self.hash_id

    def __eq__(self, other):
        return isinstance(other, Seq_Ref) and self.hash() == other.hash()

    def hash(self):
        return self.__hash__()

    def ref_dict(self):
        return dict((a, self.__getattribute__(a)) for a in ref_fieldnames)

    def setAttributes(self, version):
        title = self.title
        self.title = title if re.search('\w', title) else ''
        self.version = version
        return self
