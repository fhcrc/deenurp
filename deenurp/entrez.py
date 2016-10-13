"""
Entrez wrapper to search and fetch unlimited accessions
"""

import logging
import re
import retrying

from Bio import Entrez, SeqIO
from Bio._py3k import HTTPError
from cStringIO import StringIO

log = logging.getLogger(__name__)

RETMAX = 10000


def set_email(email):
    Entrez.email = email


def liststr(l):
    """
    Concise string representation of list.
    """
    l = list(l)
    s = ''
    if len(l) == 0:
        s = 'None'
    elif len(l) == 1:
        s = str(l[0])
    else:
        s = '{}...{}'.format(l[0], l[-1])
    return s


def entrez_pprint(prog, *extras, **kwds):
    """
    Output Entrez commands in command line format:
    http://www.ncbi.nlm.nih.gov/books/NBK179288/
    """

    entrez = dict(db='-db', rettype='-format', retmode='-mode',
                  strand='-strand', complexity='-complexity',
                  seq_start='-seq_start', seq_stop='-seq_stop', term='-term')
    keys = set(entrez.keys()).intersection(kwds.keys())
    info = [prog]
    info += list(extras)
    info += ['{} {}'.format(entrez[k], kwds[k]) for k in keys]
    return ' '.join(info)


def parse(handle, format, *args, **kwds):
    if format in ('acc', 'ft'):
        return [handle.read()]
    else:
        # SeqIO.parse does not support gbwithparts
        format = 'gb' if format == 'gbwithparts' else format
        return SeqIO.parse(handle, format, *args, **kwds)


def esearch(term, **args):
    """
    Retrieve accession numbers via ncbi search term.

    * max_records - stop retrieving records once number exceeds this value
    """

    log.info(entrez_pprint('esearch', **args))

    retstart = 0

    while True:
        try:
            handle = Entrez.esearch(term=term, retstart=retstart, **args)
            idlist = Entrez.read(handle)['IdList']

            if idlist:
                for i in idlist:
                    yield i

                retstart += args['retmax']
            else:
                return

        except HTTPError as err:
            # no need to crash everything if no results
            log.error(err)


def efetch(ids, retry=0, max_retry=10, **args):
    """
    Search GenBank given a list of accessions or GI's.
    ncbi.efetch can work unexpectedly when querying a large number of tax ids.
    The best way to mitigate is to send only its documented max of 10k ids
    at a time.  Ncbi claims it can handle more than 10k ids but will
    occassionally and unexpectedly return an http 503 error when doing so.
    """

    def print_retry_message(exception):
        """
        http exceptions: http://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html
        retrying api: https://pypi.python.org/pypi/retrying
        """

        seconds = float(retry) / 1000
        msg = '{}, retrying in {} seconds... {} max retry(ies)'.format(
            exception, seconds, max_retry)
        log.error(msg)
        return True

    @retrying.retry(
        retry_on_exception=print_retry_message,
        wait_fixed=retry,
        stop_max_attempt_number=max_retry)
    def get_records(ids=[]):
        log.info(entrez_pprint('efetch', '-id', liststr(ids), **args))
        records = parse(Entrez.efetch(id=ids, **args), args['rettype'])
        log.info('received {} {}'.format(args['rettype'], liststr(ids)))
        return records

    head, tail = 0, args['retmax']

    records = []

    while ids[head:tail]:
        records.extend(get_records(ids[head:tail]))
        head = tail
        tail += args['retmax']

    return records


def _filter_bad_records(records, ids, **args):
    passed = []
    for i, r in enumerate(records):
        print(r)
        try:
            r.decode('utf-8')
        except UnicodeError:
            log.error(entrez_pprint('efetch', '-id', ids[i], **args))
            continue
        passed.append(r)
    print(passed)
    return passed


def filter_features(records, features, strand):
    """
    http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html

    Parsing a five column, tab delimited file with a fasta file
    header line.

    example -
    > accession
    line 1
    column 1: seq_start
    column 2: seq_stop
    column 3: feature_key
    line 2
    column 4: qualifier_key
    column 5: qualifier_value
    """

    reverse_strand = {'1': '2', '2': '1'}

    # parse features for columns 3-5
    features = [f.split(':') for f in features]
    features = zip(*features)
    features = [[x if x else '.' for x in f] for f in features]
    features = ['|'.join(f) for f in features]

    # create column patterns
    column1 = '\D?(?P<seq_start>\d+)'
    column2 = '\D?(?P<seq_stop>\d+)'
    column3, column4, column5 = features
    column3 = '.*?(?P<feature_key>{})+.*?'.format(column3)
    column4 = '.*?(?P<qualifier_key>{})+.*?'.format(column4)
    column5 = '.*?(?P<qualifier_value>{})+.*?'.format(column5)

    # three line types, lines 1 and 2 are tab delimited
    accession_line = re.compile('^>.*\|(?P<accession>.*)\|.*',
                                re.IGNORECASE)
    line1 = re.compile('^{}\t{}\t{}'.format(column1, column2, column3),
                       re.IGNORECASE)
    line2 = re.compile('^\t\t\t{}\t{}'.format(column4, column5),
                       re.IGNORECASE)

    # iterate and match lines
    accession, seq_start, seq_stop = None, None, None
    for line in StringIO(records):
        match = re.search(accession_line, line)
        if match:
            accession = match.group('accession')
            continue

        match = re.search(line2, line)
        if match and seq_start and seq_stop:
            if seq_stop < seq_start:
                # swap coordinates and switch strands
                yield accession, seq_stop, seq_start, reverse_strand[strand]
            else:
                yield accession, seq_start, seq_stop, strand
            seq_start, seq_stop = None, None
            continue

        match = re.search(line1, line)
        if match:
            seq_start = int(match.group('seq_start'))
            seq_stop = int(match.group('seq_stop'))
            continue
        else:
            seq_start, seq_stop = None, None
            continue


def ffetch(features, ids, **args):
    """
    Return only specified features per list of ids
    """

    ft_args = ['db', 'retry', 'retmax', 'strand', 'retmode']
    ft_args = {k: args[k] for k in ft_args}
    ft = efetch(ids, rettype='ft', **ft_args)
    ft = filter_features(''.join(ft), features, args['strand'])
    del args['strand']
    records = []
    for accession, seq_start, seq_stop, strand in ft:
        fetched = efetch([accession], seq_start=seq_start,
                         seq_stop=seq_stop, strand=strand, **args)
        records.append((fetched[0], seq_start, seq_stop))

        # remove ids with features
        if accession in ids:
            ids.remove(accession)

    return records, ids  # ids = no features


def gbfullfetch(ids, **args):
    """
    Simple wrapper to utilize efetch while returning
    None for seq_start and seq_stop
    """

    return [(record, None, None) for record in efetch(ids, **args)], []
