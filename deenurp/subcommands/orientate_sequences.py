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
    ins = parser.add_argument_group(title='inputs')
    ins.add_argument('qseqs', help='input sequences')
    ins.add_argument('tseqs', help='target sequences')

    opt = parser.add_argument_group(title='optional inputs')
    opt.add_argument('--seq_info',
                     metavar='csv',
                     help='will split the seq_info file if provided')

    parser.add_argument('--threads',
                        metavar='NUM',
                        type=int,
                        help='number of available threads [all]')
    parser.add_argument('--id',
                        default=0.80,
                        type=float,
                        help='alignment identity percent [%(default)s]')

    # outputs
    outs = parser.add_argument_group(title='outputs')
    outs.add_argument('--out',
                      metavar='fasta',
                      default=sys.stdout,
                      help='[stdout]')
    outs.add_argument('--out_csv',
                      metavar='csv',
                      help=('output csv with columns '
                            'query,qstrand,target,id,tilo,tihi'))
    outs.add_argument('--out_notmatched',
                      metavar='fasta',
                      help='seqnames that did not match tseqs at id threshold')
    outs.add_argument('--out_notmatched_info',
                      help='seq_info file of notmatched seqs')
    outs.add_argument('--out_seq_info',
                      help='seq_info file with notmatched seqs removed')

    return parser


def read_seq_info(seq_info):
    info_dtypes = {'seqname': str, 'tax_id': str, 'tax_name': str,
                   'parent_id': str, 'rank': str, 'new_node': bool,
                   'accession': str, 'ambig_count': int}
    seq_info_df = pandas.read_csv(
        seq_info, dtype=info_dtypes,
        true_values=['yes'], false_values=['no'],
        index_col='seqname')
    return seq_info_df


def action(args):
    columns = ['query', 'target', 'qstrand', 'id', 'tilo', 'tihi']
    dtype = dict(zip(columns, [str, str, str, float, int, int]))

    prog = ['vsearch', '--usearch_global', args.qseqs,
            '--db', args.tseqs,
            '--userout', '/dev/stdout',
            '--strand', 'both',
            '--id', str(args.id),
            '--userfields', '+'.join(columns),
            '--top_hits_only',
            '--quiet']

    if args.log is not None and args.log is not sys.stdout:
        prog.extend(['--log', args.log.name])

    if args.threads:
        prog.extend(['--threads', str(args.threads)])

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
    qstrands = vsearch['qstrand'].to_dict()

    records = util.Counter(SeqIO.parse(args.qseqs, format='fasta'))
    records = (record for record in records if record.name in vsearch.index)

    with open(args.out, 'w') as out:
        for record in records:
            if qstrands[record.name] == '-':
                log.info('reversing sequence {}'.format(record.name))
                record.seq = record.seq.reverse_complement()
            SeqIO.write([record], out, 'fasta')

    if args.out_csv:
        out_columns = ['qstrand', 'target', 'id', 'tilo', 'tihi']
        vsearch.to_csv(args.out_csv, columns=out_columns)

    if args.seq_info and (args.out_seq_info or args.out_notmatched_info):
        seq_info = read_seq_info(args.seq_info)
        matched = seq_info.index.isin(vsearch.index)
        if args.out_seq_info:
            seq_info[matched].to_csv(args.out_seq_info)

        if args.out_notmatched_info:
            seq_info[~matched].to_csv(args.out_notmatched_info)
