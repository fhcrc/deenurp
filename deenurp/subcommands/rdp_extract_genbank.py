"""
Given a compressed file of GenBank records, generates a CSV file mapping from
sequence ID to tax_id, and a FASTA file with the sequences.
"""
import argparse
import csv
import logging
import gzip

from Bio import SeqIO
import sqlalchemy
from sqlalchemy.sql import select, and_
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

from deenurp.util import Counter, memoize

def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record, returning None if no taxonomy is
    available.
    """
    # Check for bad name
    try:
        source = next(i for i in gb.features if i.type == 'source')
        if 'uncultured bacterium' in ''.join(source.qualifiers.get('organism', [])).lower():
            return None
        taxon = next(i[6:] for i in source.qualifiers.get('db_xref', [])
                     if i.startswith('taxon:'))
        return taxon
    except StopIteration:
        return None

def is_classified_fn(taxonomy):
    """
    Creates a function which classifies tax_ids as classified or unclassified,
    based on presence in taxonomy and names.is_classified.
    """
    nodes = taxonomy.nodes
    @memoize
    def fetch_tax_id(tax_id):
        s = select([nodes.c.tax_id, nodes.c.is_valid, nodes.c.parent_id, nodes.c.rank])\
                .where(nodes.c.tax_id == tax_id)
        res = s.execute().fetchone()
        return res

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

def transform_id(rec):
    old_id = rec.id
    rec.id = rec.name
    rec.description = ' '.join((old_id, rec.description))
    return rec

def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)

def is_type(record):
    """
    Returns a boolean indicating whether a sequence is a member of a type strain,
    as indicated by the presence of the string '(T)' within the record description.
    """
    return '(T)' in record.description

def build_parser(p):
    p.add_argument('infile', help="""Input file, gzipped""")
    p.add_argument('database', help="""Path to taxonomy database""")
    p.add_argument('fasta_out', type=argparse.FileType('w'), help="""Path to
            write sequences in FASTA format""")
    p.add_argument('output', metavar='tax_out', type=argparse.FileType('w'),
            help="""Output path to write taxonomic information in CSV
            format""")
    p.add_argument('--no-header', action='store_false', dest='header',
            default=True, help="""Don't write a header""")

def action(a):
    """
    Run
    """
    # Start database
    e = sqlalchemy.create_engine('sqlite:///{0}'.format(a.database))
    tax = Taxonomy(e, ncbi.ranks)
    is_classified = is_classified_fn(tax)

    with gzip.open(a.infile) as fp, a.output as out_fp, a.fasta_out as fasta_fp:
        records = Counter(SeqIO.parse(fp, 'genbank'), prefix='Record ')
        taxa = ((record, tax_of_genbank(record)) for record in records)
        taxa = ((i, j if not j or is_classified(j) else None) for i, j in taxa)

        writer = csv.writer(out_fp, lineterminator='\n',
                quoting=csv.QUOTE_MINIMAL)
        if a.header:
            writer.writerow(('seqname', 'tax_id', 'accession', 'description',
                'length', 'ambig_count', 'is_type', 'rdp_lineage'))
        for record, tax_id in taxa :
            accession = record.id
            row = (record.name, tax_id, accession, record.description,
                    len(record), count_ambiguous(str(record.seq)),
                    str(is_type(record)).upper(),
                    ';'.join(record.annotations.get('taxonomy', [])).replace('"', ''))
            writer.writerow(row)
            SeqIO.write([transform_id(record)], fasta_fp, 'fasta')

    logging.info("Total records: %d", records.count)
