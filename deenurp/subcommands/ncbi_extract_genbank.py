#!/usr/bin/env python

#    Entrez.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Entrez.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Entrez.py.  If not, see <http://www.gnu.org/licenses/>.

"""
NCBI Entrez tool for downloading nucleotide Genbank or Accession numbers
given list of accession numbers, GI numbers or esearch -term.  Provides
multiprocessing and record chunking.

general rules: http://www.ncbi.nlm.nih.gov/books/NBK25497/
esearch and efetch: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4
sequence identifiers: http://www.ncbi.nlm.nih.gov/genbank/sequenceids/
feature tables: http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html
http retrying: https://pypi.python.org/pypi/retrying
"""

import argparse
import functools
import itertools
import logging
import multiprocessing
import pandas
import re
import sys

from deenurp import util, entrez
from sqlalchemy.sql import select

log = logging.getLogger(__name__)

seq_info_columns = ['seqname', 'version', 'accession', 'name', 'description',
                    'gi', 'tax_id', 'date', 'source', 'keywords', 'organism',
                    'length', 'ambig_count', 'is_type', 'taxid_classified']

ref_info_columns = ['version', 'title', 'authors', 'comment',
                    'consrtm', 'journal', 'medline_id', 'pubmed_id']


def build_parser(parser):
    parser.add_argument('email', help='user email')

    parser.add_argument('ids',
                        type=util.file_opener('r'),
                        help=('file list of GIs or accessions to fetch; '
                              'will be added to any search terms'))

    parser.add_argument('--taxonomy',
                        metavar='DB',
                        help=('Path to taxonomy database for '
                              'validating and updating tax_ids'))
    parser.add_argument('--threads',
                        metavar='NUM',
                        default=int(multiprocessing.cpu_count()),
                        type=int,
                        help='number of available threads [%(default)s]')
    parser.add_argument('--chunksize',
                        metavar='NUM',
                        default=entrez.RETMAX,
                        type=int,
                        help="""number of records to return
                                per query max 10k [%(default)s]""")
    parser.add_argument('--retry',
                        metavar='MILLISECONDS',
                        type=int,
                        default=60000,
                        help="""after http exception time to
                                wait before trying again [%(default)s]""")
    parser.add_argument('--max-records',
                        type=int,
                        metavar='N',
                        help="""limit number of records from
                                esearch or ids file to N""")

    parser.add_argument('--all-versions',
                        action='store_true',
                        help='retrieve all given record versions')
    parser.add_argument('--strand',
                        default='1',
                        choices=['1', '2'],
                        help="""1 for the plus strand, 2 for
                                the minus strand. [%(default)s]""")
    parser.add_argument('--feature',
                        action='append',
                        dest='features',
                        help="""parse only specific record features in feature
                                table columns via string pattern
                                feature_key:qualifier_key:qualifier_value.
                                ex: rRNA:product:16S""")

    outs = parser.add_argument_group('various output options')
    outs.add_argument('--out',
                      metavar='FASTA',
                      type=argparse.FileType('w'),
                      default=sys.stdout,
                      help='[stdout]')
    outs.add_argument('--info-out',
                      metavar='CSV',
                      help='write seq_info file')
    outs.add_argument('--refs-out',
                      metavar='CSV',
                      help='output references')

    return parser


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
def fetch_tax_info(tax_id, taxonomy):
    c = taxonomy.nodes.c  # columns
    s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
    s = s.where(c.tax_id == tax_id)
    return s.execute().fetchone()


@util.memoize
def species_is_classified(tax_id, taxonomy):
    """
    return is_valid from ncbi taxonomy
    """
    tax_id, is_valid, parent_id, rank = fetch_tax_info(tax_id, taxonomy)
    if rank == 'species':
        return is_valid
    elif tax_id == parent_id:
        return False
    else:
        # instead of recursion create a table of ranks by index. code
        # will be something liek ranks.index(rank) > species_index
        return species_is_classified(parent_id, taxonomy)


@util.memoize
def update_taxid(tax_id, name, taxonomy):
    """
    Update the tax_id if it has been merged.  If the tax_id is not in the
    database a tax_id will be assigned based on the organism name
    """
    try:
        taxonomy._node(tax_id)
    except KeyError as err:
        new_tax_id = taxonomy._get_merged(tax_id)
        if new_tax_id != tax_id:
            msg = 'tax_id {} merged to {}'.format(tax_id, new_tax_id)
            logging.warn(msg)
            tax_id = new_tax_id
        elif name:
            try:
                tax_id, _, _ = taxonomy.primary_from_name(name)
            except KeyError as err:
                logging.warn(err)
                tax_id = None
        else:
            msg = 'taxid {} not found in taxonomy, dropping'
            logging.warn(msg.format(tax_id))
            tax_id = None

    return tax_id


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def latest_versions(ids):
    """
    http://www.ncbi.nlm.nih.gov/genbank/sequenceids/
    """

    ordered = sorted(ids, reverse=True)
    versions = lambda x: re.search('.+\.?', x).group()
    grouped = itertools.groupby(ordered, key=versions)

    for _, ids in grouped:
        yield ids.next()  # top id is latest record version


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


def parse_record(record, taxonomy=None, seq_start=None,
                 seq_stop=None, strand=None):
    accession, version = accession_version_of_genbank(record)
    organism = record.annotations['organism']

    tax_id = tax_of_genbank(record)
    if taxonomy is not None:
        tax_id = update_taxid(tax_id, organism, taxonomy)
        is_classified = species_is_classified(tax_id, taxonomy)
    else:
        is_classified = False

    seq = str(record.seq)

    if None not in (seq_start, seq_stop):
        seq_id = '{}_{}_{}'.format(accession, seq_start, seq_stop)
    else:
        seq_id = record.id

    info = dict(accession=accession,
                ambig_count=count_ambiguous(seq),
                length=len(record),
                seqname=seq_id,
                is_type=is_type(record),
                date=record.annotations['date'],
                description=record.description,
                keywords=';'.join(record.annotations.get('keywords', [])),
                name=record.name,
                organism=organism,
                seq=seq,
                source=record.annotations['source'],
                tax_id=tax_id,
                taxid_classified=is_classified,
                seq_start=seq_start,
                seq_stop=seq_stop,
                strand=strand,
                version=version)
    return pandas.Series(info)


def parse_references(record, version):
    references = []
    for r in record.annotations['references']:
        references.append(
            dict(version=version,
                 title=r.title,
                 authors=r.authors,
                 comment=r.comment,
                 consrtm=r.consrtm,
                 journal=r.journal,
                 medline_id=r.medline_id,
                 pubmed_id=r.pubmed_id))
    return pandas.DataFrame(references)


def action(args):
    entrez.set_email(args.email)

    chunksize = min(entrez.RETMAX, args.chunksize)  # 10k is ncbi max right now

    ids = (i.strip() for i in args.ids)  # remove newlines, etc
    ids = itertools.chain(ids, args.ids)
    ids = itertools.islice(ids, args.max_records)

    # take latest version accession
    if not args.all_versions:
        log.info('returning only latest versions of ids')
        ids = latest_versions(ids)

    # break up ids into specified chunks
    ids = (chunk for chunk in util.chunker(ids, chunksize))

    fetch_args = dict(db='nucleotide',
                      retry=args.retry,
                      rettype='gbwithparts',
                      retmax=chunksize,
                      retmode='text',
                      strand=args.strand,
                      complexity=1)

    if args.features:
        func = functools.partial(entrez.ffetch, args.features, **fetch_args)
    else:
        func = functools.partial(entrez.gbfullfetch, **fetch_args)

    records, references = [], []
    pool = multiprocessing.Pool(processes=args.threads)
    try:
        gbs = itertools.chain.from_iterable(pool.imap_unordered(func, ids))
    except Exception as e:
        log.error(e)

    for g, seq_start, seq_stop in gbs:
        record = parse_record(g, taxonomy=args.taxonomy,
                              seq_start=seq_start, seq_stop=seq_stop)
        records.append(record)
        if args.refs_out:
            references.append(parse_references(g, record['version']))

    df_records = pandas.DataFrame(records)
    to_fasta = lambda x: '>{seqname}\n{seq}'.format(**x)
    args.out.write('\n'.join(df_records.apply(to_fasta, axis=1)))

    if args.info_out:
        df_records.to_csv(
            args.info_out, index=False, columns=seq_info_columns)

    if args.refs_out:
        pandas.concat(references).to_csv(
            args.refs_out, index=False, columns=ref_info_columns)
