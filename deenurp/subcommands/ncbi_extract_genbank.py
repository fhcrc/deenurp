"""
NCBI Entrez tool for downloading nucleotide Genbank records
given list of accession numbers, GI numbers. Provides
multiprocessing and record chunking.

general rules: http://www.ncbi.nlm.nih.gov/books/NBK25497/
esearch and efetch: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4
sequence identifiers: http://www.ncbi.nlm.nih.gov/genbank/sequenceids/
feature tables: http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html
http retrying: https://pypi.python.org/pypi/retrying
"""

import csv
import functools
import itertools
import logging
import multiprocessing
import re
import sqlalchemy

from deenurp import util, entrez
from sqlalchemy.sql import select
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

log = logging.getLogger(__name__)

seq_info_columns = ['seqname', 'version', 'accession', 'name', 'description',
                    'gi', 'tax_id', 'date', 'source', 'keywords', 'organism',
                    'length', 'ambig_count', 'is_type', 'taxid_classified',
                    'seq_start', 'seq_stop']

ref_info_columns = ['version', 'title', 'authors', 'comment',
                    'consrtm', 'journal', 'medline_id', 'pubmed_id']


def build_parser(parser):
    # ins
    parser.add_argument('email', help='user email')
    parser.add_argument('ids',
                        type=util.file_opener(mode='r'),
                        help=('file list of GIs or accessions to fetch; '
                              'will be added to any search terms'))

    # outs
    parser.add_argument('fasta_out',
                        metavar='FASTA',
                        type=util.file_opener(mode='w'),
                        help='[stdout]')
    parser.add_argument('info_out',
                        metavar='CSV',
                        type=util.file_opener(mode='w'),
                        help='write seq_info file')
    parser.add_argument('--refs-out',
                        metavar='CSV',
                        type=util.file_opener(mode='w'),
                        help='output references')

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
            log.warn(msg)
            tax_id = new_tax_id
        elif name:
            try:
                tax_id, _, _ = taxonomy.primary_from_name(name)
            except KeyError as err:
                log.warn(err)
                tax_id = None
        else:
            msg = 'taxid {} not found in taxonomy, dropping'
            log.warn(msg.format(tax_id))
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
                gi=record.annotations.get('gi', ''),
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
    return info


def parse_references(record):
    references = []
    for r in record.annotations['references']:
        references.append(
            dict(title=r.title,
                 authors=r.authors,
                 comment=r.comment,
                 consrtm=r.consrtm,
                 journal=r.journal,
                 medline_id=r.medline_id,
                 pubmed_id=r.pubmed_id))
    return references


def action(args):
    entrez.set_email(args.email)

    chunksize = min(entrez.RETMAX, args.chunksize)  # 10k is ncbi max right now

    ids = (i.strip() for i in args.ids)  # remove newlines, etc
    ids = itertools.islice(ids, args.max_records)

    # take latest version accession
    if not args.all_versions:
        log.info('returning only latest versions of ids')
        ids = latest_versions(ids)

    # break up ids into specified chunks
    ids = (chunk for chunk in util.chunker(ids, chunksize))

    if args.taxonomy:
        e = sqlalchemy.create_engine('sqlite:///{0}'.format(args.taxonomy))
        taxonomy = Taxonomy(e, ncbi.ranks)
    else:
        taxonomy = None

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

    pool = multiprocessing.Pool(processes=args.threads)
    gbs = itertools.chain.from_iterable(pool.imap_unordered(func, ids))
    # gbs = itertools.chain.from_iterable(itertools.imap(func, ids))  # testing

    info_out = csv.DictWriter(
        args.info_out, fieldnames=seq_info_columns, extrasaction='ignore')
    info_out.writeheader()
    if args.refs_out:
        refs_out = csv.DictWriter(args.refs_out, fieldnames=ref_info_columns)
        refs_out.writeheader()

    for g, seq_start, seq_stop in util.Counter(gbs, report_every=0.01):
        record = parse_record(
            g, taxonomy=taxonomy, seq_start=seq_start, seq_stop=seq_stop)
        args.fasta_out.write('>{seqname}\n{seq}\n'.format(**record))
        info_out.writerow(record)
        if args.refs_out:
            # add version
            version = record['version']
            refs = [dict(version=version, **r) for r in parse_references(g)]
            refs_out.writerows(refs)
