"""
Fix orientation of sequences and output target sequence alignment indexes
"""

import logging
import pandas
import subprocess
import sys

from Bio import SeqIO
from cStringIO import StringIO
from deenurp import util

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('qseqs', metavar='fasta', help='input sequences')
    parser.add_argument('tseqs', metavar='fasta', help='target sequences')

    parser.add_argument('--threads',
                        metavar='NUM',
                        type=int,
                        help='number of available threads [all]')
    parser.add_argument('--id',
                        default=.8,
                        type=float,
                        help='alignment identity percent')

    # outputs
    outs = parser.add_argument_group(title='outputs')
    outs.add_argument('--out',
                      metavar='fasta',
                      default=sys.stdout,
                      help='[stdout]')
    outs.add_argument('--out_csv',
                      metavar='csv',
                      help='output csv with columns query,target,tilo,tihi')
    outs.add_argument('--out_notmatched',
                      metavar='fasta',
                      help='seqnames that did not match tseqs at id threshold')
    return parser


def action(args):
    columns = ['query', 'target', 'qstrand', 'tilo', 'tihi']
    dtype = dict(zip(columns, [str, str, str, int, int]))

    prog = ['vsearch', '--usearch_global', args.qseqs,
            '--db', args.tseqs,
            '--userout', '/dev/stdout',
            '--strand', 'both',
            '--id', str(args.id),
            '--userfields', '+'.join(columns),
            '--top_hits_only',
            '--quiet']

    if args.threads:
        prog.extend(['--threads', args.threads])

    if args.out_notmatched:
        prog.extend(['--notmatched', args.out_notmatched])

    log.info(' '.join(prog))

    pipe = subprocess.Popen(
        prog, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    results, errors = pipe.communicate()

    errors = errors.strip()
    if errors:
        log.error(errors)

    vsearch = pandas.read_csv(
        StringIO(results), dtype=dtype, names=columns,
        header=None, sep='\t', index_col='query')

    records = util.Counter(SeqIO.parse(args.qseqs, format='fasta'))
    records = (record for record in records if record.name in vsearch.index)
    with open(args.out, 'w') as out:
        for record in records:
            row = vsearch[vsearch.index == record.name]
            if row[row['qstrand'] == '+'].empty:
                log.info('reversing sequence {}'.format(record.name))
                record.seq = record.seq.reverse_complement()
            SeqIO.write([record], out, 'fasta')

    if args.out_csv:
        vsearch.to_csv(args.out_csv, columns=['target', 'tilo', 'tihi'])
