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
NCBI Entrez tool for downloading nucleotide Accession numbers
given esearch -term.  Provides multiprocessing and record chunking.

general rules: http://www.ncbi.nlm.nih.gov/books/NBK25497/
esearch and efetch: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4
sequence identifiers: http://www.ncbi.nlm.nih.gov/genbank/sequenceids/
http retrying: https://pypi.python.org/pypi/retrying
"""

import argparse
import functools
import itertools
import logging
import multiprocessing
import sys

from deenurp import entrez, util

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('email', help='user email')
    parser.add_argument('term', help='esearch term')

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
                        help=('limit number of ids from esearch term'))
    parser.add_argument('--all-versions',
                        action='store_true',
                        help='retrieve all given record versions')

    # outputs
    parser.add_argument('--out',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='[stdout]')

    return parser


def action(args):
    entrez.set_email(args.email)
    chunksize = min(entrez.RETMAX, args.chunksize)  # 10k is ncbi max right now
    search_args = dict(retmax=chunksize, db='nucleotide', rettype='uilist')
    ids = entrez.esearch(args.term, **search_args)
    ids = itertools.islice(ids, args.max_records)
    id_chunks = (chunk for chunk in util.chunker(ids, chunksize))
    fetch_args = dict(db='nucleotide', retry=args.retry, rettype='acc',
                      retmax=chunksize, retmode='text', complexity=1)
    func = functools.partial(entrez.efetch, **fetch_args)
    pool = multiprocessing.Pool(processes=args.threads)
    for versions in pool.imap_unordered(func, id_chunks):
        # chunking on RETMAX so we can just write the top result
        args.out.write(versions[0])
