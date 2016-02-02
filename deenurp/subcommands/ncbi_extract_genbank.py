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
import os
import re

from deenurp import util, entrez

log = logging.getLogger(__name__)

seq_info_columns = ['seqname', 'version', 'accession', 'name',
                    'description', 'gi', 'tax_id', 'date', 'source',
                    'keywords', 'organism', 'length', 'ambig_count',
                    'is_type', 'seq_start', 'seq_stop']

reference_columns = ['pubmed_id', 'medline_id', 'title',
                     'authors', 'journal', 'consrtm', 'comment']

pubmed_columns = ['pubmed_id', 'version', 'accession']


def build_parser(parser):
    # ins
    parser.add_argument('email', help='user email')
    parser.add_argument('ids',
                        type=util.file_opener(mode='r'),
                        help=('file list of GIs or accessions to fetch; '
                              'will be added to any search terms'))

    # outs
    parser.add_argument('fasta_out',
                        metavar='fasta',
                        help='sequence file')
    parser.add_argument('info_out',
                        metavar='csv',
                        help='write seq_info file')
    parser.add_argument('--pubmed_ids',
                        metavar='csv',
                        help=('csv with columns '
                              '[version, accession, pubmed_id]'))
    parser.add_argument('--references',
                        metavar='csv',
                        help=('reference details'))
    parser.add_argument('--no-features',
                        type=util.file_opener(mode='a'),
                        help=('output version numbers of records with no '
                              'matching features'))

    parser.add_argument('--threads',
                        metavar='int',
                        default=int(multiprocessing.cpu_count()),
                        type=int,
                        help='number of available threads [%(default)s]')
    parser.add_argument('--chunksize',
                        metavar='int',
                        default=entrez.RETMAX,
                        type=int,
                        help="""number of records to return
                                per query max 10k [%(default)s]""")
    parser.add_argument('--retry',
                        metavar='milliseconds',
                        type=int,
                        default=60000,
                        help="""after http exception time to
                                wait before trying again [%(default)s]""")
    parser.add_argument('--max-records',
                        type=int,
                        metavar='int',
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
    parser.add_argument('--type-strains',
                        type=util.file_opener(mode='r'),
                        help='list of type strain records by gi numbers')

    return parser


def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record,
    returning None if no taxonomy is present
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


def parse_record(record, seq_start=None, seq_stop=None,
                 strand=None, types=set()):
    accession, version = accession_version_of_genbank(record)
    organism = record.annotations['organism']

    tax_id = tax_of_genbank(record)
    seq = str(record.seq)

    if None not in (seq_start, seq_stop):
        seq_id = '{}_{}_{}'.format(accession, seq_start, seq_stop)
    else:
        seq_id = record.id

    info = dict(accession=accession,
                ambig_count=count_ambiguous(seq),
                length=len(record),
                seqname=seq_id,
                is_type=version in types,
                date=record.annotations['date'],
                description=record.description,
                gi=record.annotations.get('gi', ''),
                keywords=';'.join(record.annotations.get('keywords', [])),
                name=record.name,
                organism=organism,
                seq=seq,
                source=record.annotations['source'],
                tax_id=tax_id,
                seq_start=seq_start,
                seq_stop=seq_stop,
                strand=strand,
                version=version)
    return info


def parse_references(record):
    """
    Parse reference annotations that have a pubmed_id
    """
    references = []
    if 'references' in record.annotations:
        refs = [r for r in record.annotations['references'] if r.pubmed_id]
        for r in refs:
            references.append(
                dict(title=r.title,
                     authors=r.authors,
                     comment=r.comment,
                     consrtm=r.consrtm,
                     journal=r.journal,
                     medline_id=r.medline_id,
                     pubmed_id=r.pubmed_id))
    return references


def write_or_append(fasta_path, info_path, pubmed_ids_path, references_path):
    if os.path.isfile(fasta_path):
        log.info('{} exists, appending new results'.format(fasta_path))
        fasta_out = open(fasta_path, 'a')
    else:
        fasta_out = open(fasta_path, 'w')

    if os.path.isfile(info_path):
        log.info('{} exists, appending new results'.format(info_path))

        # get column order
        with open(info_path) as info_in:
            fieldnames = csv.DictReader(info_in).fieldnames
        if set(fieldnames) != set(seq_info_columns):
            msg = 'input columns {} != {}'.format(fieldnames, seq_info_columns)
            raise ValueError(msg)

        info_out = csv.DictWriter(
            open(info_path, 'a'), fieldnames=fieldnames, extrasaction='ignore')
    else:
        info_out = csv.DictWriter(
            open(info_path, 'w'),
            fieldnames=seq_info_columns,
            extrasaction='ignore')
        info_out.writeheader()

    if pubmed_ids_path is not None and os.path.isfile(pubmed_ids_path):
        log.info('{} exists, appending new results'.format(pubmed_ids_path))

        # get column order
        with open(pubmed_ids_path) as pubmed_ids_in:
            fieldnames = csv.DictReader(pubmed_ids_in).fieldnames
        if set(fieldnames) != set(pubmed_columns):
            msg = 'input columns {} != {}'.format(fieldnames, pubmed_columns)
            raise ValueError(msg)

        pubmed_ids_out = csv.DictWriter(
            open(pubmed_ids_path, 'a'), fieldnames=fieldnames)
    elif pubmed_ids_path is not None:
        pubmed_ids_out = csv.DictWriter(
            open(pubmed_ids_path, 'w'), fieldnames=pubmed_columns)
        pubmed_ids_out.writeheader()
    else:
        pubmed_ids_out = None

    if references_path is not None and os.path.isfile(references_path):
        log.info('{} exists, appending new results'.format(references_path))

        # get column order
        with open(references_path) as references_in:
            fieldnames = csv.DictReader(references_in).fieldnames
        if set(fieldnames) != set(reference_columns):
            msg = 'input columns {} != {}'
            raise ValueError(msg.format(fieldnames, reference_columns))

        references_out = csv.DictWriter(
            open(references_path, 'a'), fieldnames=fieldnames)
    elif references_path is not None:
        references_out = csv.DictWriter(
            open(references_path, 'w'), fieldnames=reference_columns)
        references_out.writeheader()
    else:
        references_out = None

    return fasta_out, info_out, pubmed_ids_out, references_out


def parse_type_strain_set(types_file):
    if types_file:
        types = (t.strip() for t in types_file)
        types = set(t for t in types if t)
    else:
        types = set()

    return types


def expand_accession_ranges(accessions):
    for a in accessions:
        if '-' in a or ':' in a:
            prefix = re.findall('\D+', a)[0]
            a_range = re.findall('\d+', a)
            a_range = map(int, a_range)  # integers for xrange
            a_range = xrange(*a_range)  # expand range
            a_range = map(str, a_range)  # back to strings
            for acc in a_range:
                yield prefix + acc
        else:
            yield a


def action(args):
    entrez.set_email(args.email)

    chunksize = min(entrez.RETMAX, args.chunksize)  # 10k is ncbi max right now

    ids = (i for i in args.ids if not i.startswith('#'))  # ignore comments
    ids = (i.strip() for i in ids)  # remove newlines
    ids = (i for i in ids if i)  # ignore blank lines
    ids = expand_accession_ranges(ids)  # id ranges
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

    pool = multiprocessing.Pool(processes=args.threads)
    records = pool.imap_unordered(func, ids)
    # records = itertools.imap(func, ids)

    fasta_out, info_out, pubmed_ids_out, references_out = write_or_append(
        args.fasta_out, args.info_out, args.pubmed_ids, args.references)

    types = parse_type_strain_set(args.type_strains)

    for gbs, no_features in records:
        for g, seq_start, seq_stop in gbs:
            record = parse_record(
                g, seq_start=seq_start, seq_stop=seq_stop, types=types)
            fasta_out.write('>{seqname}\n{seq}\n'.format(**record))
            info_out.writerow(record)

            if pubmed_ids_out:
                references = parse_references(g)
                for r in references:
                    pubmed_ids_out.writerow(
                        {'version': record['version'],
                         'accession': record['accession'],
                         'pubmed_id': r['pubmed_id']})
                if references_out:
                    references_out.writerows(references)

        if no_features:
            log.warn('no features found ' + entrez.liststr(no_features))
            if args.no_features:
                args.no_features.write('\n'.join(no_features) + '\n')
