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

from deenurp import util, seq_info

type_keywords = ['(T)', 'ATCC', 'NCTC', 'NBRC', 'CCUG',
                 'DSM', 'JCM', 'NCDO', 'NCIB', 'CIP']


def species_is_classified_fn(taxonomy):
    """
    Creates a function which classifies tax_ids as classified or
    unclassified at species level according to taxonomy database,
    based on presence in taxonomy and names.is_classified.
    """
    nodes = taxonomy.nodes
    ranks = taxonomy.ranks
    species_index = ranks.index('species')

    @util.memoize
    def fetch_tax_id(tax_id):
        # get tax node data
        c = nodes.c
        s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
        s = s.where(c.tax_id == tax_id)
        return s.execute().fetchone()

    @util.memoize
    def is_classified(tax_id):
        res = fetch_tax_id(tax_id)
        if not res:
            return False

        tax_id, is_valid, parent_id, rank = res
        if rank == 'species':
            return is_valid
        elif rank == 'no rank':
            # NOTE: is this only root?  Recusively move up taxonomy,
            # if above species will return False.
            if tax_id == parent_id:
                return False
            return is_classified(parent_id)
        else:
            # FIXME: ranks needs to be imported from ncbi_taxonomy.db database
            # returns FALSE if above species
            return ranks.index(rank) > species_index

    return is_classified


def transform_id(rec):
    old_id = rec.id
    rec.id = rec.name
    rec.description = ' '.join((old_id, rec.description))
    return rec


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def update_taxid(tax_id, taxonomy, name):
    """
    """

    @util.memoize
    def update(tax_id, name):
        try:
            taxonomy._node(tax_id)
        except KeyError as err:
            new_tax_id = taxonomy._get_merged(tax_id)
            if new_tax_id != tax_id:
                msg = 'updating tax_id {} to {}'.format(tax_id, new_tax_id)
                logging.warn(msg)
                tax_id = new_tax_id
            elif name:
                try:
                    tax_id, _, _ = taxonomy.primary_from_name(name)
                except KeyError as err:
                    logging.warn(err)
                    tax_id = None
            else:
                msg = 'taxid {} not found in taxonomy, dropping'.format(tax_id)
                logging.warn(msg)
                tax_id = None

        return tax_id

    return update(tax_id, name)


def build_parser(p):
    p.add_argument('infile', type=util.file_opener('r'),
                   help="""Input file, gzipped""")
    p.add_argument('database', help="""Path to taxonomy database""")
    p.add_argument('fasta_out', type=util.file_opener('w'),
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
    taxonomy = Taxonomy(e, ncbi.ranks)
    is_classified = species_is_classified_fn(taxonomy)

    with a.infile as fp, \
            a.output as out_fp, \
            a.fasta_out as fasta_fp:
        records = SeqIO.parse(fp, 'genbank')
        records = util.Counter(records, prefix='Record ')
        taxa = ((record, seq_info.tax_of_genbank(record))
                for record in records)
        taxa = ((record, tax_id, record.annotations['organism'])
                for record, tax_id in taxa)
        taxa = ((record, update_taxid(tax_id, taxonomy, organism))
                for record, tax_id, organism in taxa)

        writer = csv.writer(out_fp,
                            lineterminator='\n',
                            quoting=csv.QUOTE_NONNUMERIC)
        if a.header:
            header = ('version', 'seqname', 'tax_id', 'accession',
                      'description', 'length', 'ambig_count',
                      'is_type', 'rdp_lineage', 'taxid_classified')
            writer.writerow(header)

        for record, tax_id in taxa:
            accession, version = seq_info.accession_version_of_genbank(record)
            rdp_lineage = ';'.join(record.annotations.get('taxonomy', []))
            rdp_lineage = rdp_lineage.replace('"', '')

            row = (version, record.name, tax_id, accession, record.description,
                   len(record), count_ambiguous(str(record.seq)),
                   str(seq_info.is_type(record)).upper(),
                   rdp_lineage, is_classified(tax_id))

            writer.writerow(row)
            SeqIO.write([transform_id(record)], fasta_fp, 'fasta')

    logging.info("Total records: %d", records.count)
