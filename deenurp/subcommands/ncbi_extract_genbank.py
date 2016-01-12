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
                        help='')
    parser.add_argument('info_out',
                        metavar='CSV',
                        help='write seq_info file')
    parser.add_argument('refs_out',
                        metavar='CSV',
                        help='output references')
    parser.add_argument('--no-features',
                        type=util.file_opener(mode='a'),
                        help=('output version numbers of records with no '
                              'matching features'))

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
                is_type=record in types,
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
    references = []
    if 'references' in record.annotations:
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


def write_or_append(info_path, refs_path, fasta_path):
    if os.path.isfile(info_path):
        log.info('{} exists, appending new results'.format(info_path))

        # get column order
        with open(info_path) as info_out:
            fieldnames = csv.DictReader(info_out).fieldnames
        if set(fieldnames) != set(seq_info_columns):
            msg = 'input columns {} != {}'.format(fieldnames, seq_info_columns)
            raise ValueError(msg)

        info_out = csv.DictWriter(
            open(info_path, 'a'),
            fieldnames=fieldnames,
            extrasaction='ignore')
    else:
        info_out = csv.DictWriter(
            open(info_path, 'w'),
            fieldnames=seq_info_columns,
            extrasaction='ignore')
        info_out.writeheader()

    if os.path.isfile(refs_path):
        log.info('{} exists, appending new results'.format(refs_path))

        # get column order
        with open(refs_path) as refs_out:
            fieldnames = csv.DictReader(refs_out).fieldnames
        if set(fieldnames) != set(ref_info_columns):
            msg = 'input columns {} != {}'.format(fieldnames, ref_info_columns)
            raise ValueError(msg)

        refs_out = csv.DictWriter(open(refs_path, 'a'), fieldnames=fieldnames)
    else:
        refs_out = csv.DictWriter(
            open(refs_path, 'w'), fieldnames=ref_info_columns)
        refs_out.writeheader()

    if os.path.isfile(fasta_path):
        log.info('{} exists, appending new results'.format(fasta_path))
        fasta_out = open(fasta_path, 'a')
    else:
        fasta_out = open(fasta_path, 'w')

    return info_out, refs_out, fasta_out


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
    ids = expand_accession_ranges(ids)  # ranges
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

    info_out, refs_out, fasta_out = write_or_append(
        args.info_out, args.refs_out, args.fasta_out)

    types = parse_type_strain_set(args.type_strains)

    for gbs, no_features in records:
        for g, seq_start, seq_stop in gbs:
            record = parse_record(
                g, seq_start=seq_start, seq_stop=seq_stop, types=types)
            fasta_out.write('>{seqname}\n{seq}\n'.format(**record))
            info_out.writerow(record)

            version = record['version']
            refs = [dict(version=version, **r)  # add version
                    for r in parse_references(g)]
            refs_out.writerows(refs)

        if no_features:
            log.warn('no features found ' + entrez.liststr(no_features))
            if args.no_features:
                args.no_features.write('\n'.join(no_features) + '\n')
