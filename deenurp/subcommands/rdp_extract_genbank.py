"""
Extract sequences and annotation from a file of GenBank records.
"""

import argparse
import csv
import logging

from Bio import SeqIO
import sqlalchemy
from sqlalchemy.sql import select
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

from deenurp.util import (Counter, memoize, file_opener,
                          accession_version_of_genbank,
                          tax_of_genbank)

type_keywords = ['(T)', 'ATCC', 'NCTC', 'NBRC', 'CCUG',
                 'DSM', 'JCM', 'NCDO', 'NCIB', 'CIP']


def is_type(record):
    """
    Returns a boolean indicating whether a sequence is a member of a type
    strain, as indicated by the presence of the string '(T)' within the
    record description.
    """
    for t in type_keywords:
        if t in record.description:
            return True

    return False


def is_classified_fn(taxonomy):
    """
    Creates a function which classifies tax_ids as classified or unclassified,
    based on presence in taxonomy and names.is_classified.
    """
    nodes = taxonomy.nodes

    @memoize
    def fetch_tax_id(tax_id):
        # get tax node data
        c = nodes.c
        s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
        s = s.where(c.tax_id == tax_id)
        return s.execute().fetchone()

    @memoize
    def is_classified(tax_id):
        res = fetch_tax_id(tax_id)
        if not res:
            return False

        tax_id, is_class, parent_id, rank = res
        if rank == 'species':
            return is_class
        else:
            if tax_id == parent_id:
                return False
            return is_classified(parent_id)

    return is_classified


def update_taxid(tax_id, taxonomy):
    """
    Check the taxonomy for latest taxid and return it
    """
    merged = taxonomy.merged

    # fetch updated tax_id if exists
    c = merged.c
    s = select([c.new_tax_id]).where(c.old_tax_id == tax_id)
    new_tax_id = s.execute().fetchone()
    return new_tax_id[0] if new_tax_id else tax_id


def transform_id(rec):
    old_id = rec.id
    rec.id = rec.name
    rec.description = ' '.join((old_id, rec.description))
    return rec


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def build_parser(p):
    p.add_argument('infile', type=file_opener('r'),
                   help="""Input file, gzipped""")
    p.add_argument('database', help="""Path to taxonomy database""")
    p.add_argument('fasta_out', type=file_opener('w'),
                   help="""Path to write sequences in FASTA format.
                           Specify '.gz' or '.bz2' extension to compress.""")
    p.add_argument('output', metavar='tax_out', type=argparse.FileType('w'),
                   help="""Output path to write taxonomic
                           information in CSV format""")
    p.add_argument('--no-header', action='store_false', dest='header',
                   default=True, help="""Don't write a header""")


def action(a):
    # Start database
    e = sqlalchemy.create_engine('sqlite:///{0}'.format(a.database))
    tax = Taxonomy(e, ncbi.ranks)
    is_classified = is_classified_fn(tax)

    with a.infile as fp, \
            a.output as out_fp, \
            a.fasta_out as fasta_fp:
        records = SeqIO.parse(fp, 'genbank')
        records = Counter(records, prefix='Record ')
        taxa = ((rec, tax_of_genbank(rec)) for rec in records)
        taxa = ((i, update_taxid(j, tax)) for i, j in taxa)

        writer = csv.writer(out_fp,
                            lineterminator='\n',
                            quoting=csv.QUOTE_NONNUMERIC)
        if a.header:
            header = ('version', 'seqname', 'tax_id', 'accession',
                      'description', 'length', 'ambig_count',
                      'is_type', 'rdp_lineage', 'taxid_classified')
            writer.writerow(header)

        for record, tax_id in taxa:
            accession, version = accession_version_of_genbank(record)
            rdp_lineage = ';'.join(record.annotations.get('taxonomy', []))
            rdp_lineage = rdp_lineage.replace('"', '')

            row = (version, record.name, tax_id, accession, record.description,
                   len(record), count_ambiguous(str(record.seq)),
                   str(is_type(record)).upper(),
                   rdp_lineage, is_classified(tax_id))

            writer.writerow(row)
            SeqIO.write([transform_id(record)], fasta_fp, 'fasta')

    logging.info("Total records: %d", records.count)
