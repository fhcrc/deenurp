"""
deduplicate sequences using md5 sum as sequence names
"""
import csv
import hashlib
import logging
import pandas

from Bio import SeqIO

from deenurp import util

log = logging.getLogger(__name__)


def build_parser(p):
    # inputs
    p.add_argument(
        'sequences',
        help="""sequence file""")
    p.add_argument(
        'seq_info',
        help="""Sequence metadata""")

    # options
    p.add_argument(
        '--group-by',
        metavar='columns',
        help=('comma separated list of columns to '
              'group on and calculate weights'))
    p.add_argument(
        '--prefer-columns',
        metavar='columns',
        help=('columns to sort by and choose as seq_info representative'))

    # outputs
    seq_outs = p.add_argument_group('outputs for sequences')
    seq_outs.add_argument(
        'out',
        help='dedup sequences with md5 as seqname')
    seq_outs.add_argument(
        'out_info',
        help='output seqinfo file with new seqhash column')


def action(args):
    dtype = {'gi': str, 'tax_id': str, 'species': str}
    seq_info = pandas.read_csv(args.seq_info, dtype=dtype, index_col='seqname')

    log.info('reading sequences')
    with util.file_opener()(args.sequences) as sequences_in:
        seqhashes = dict()
        for record in util.Counter(SeqIO.parse(sequences_in, 'fasta')):
            seq = str(record.seq).replace('\n', '').upper()
            seqhashes[record.name] = hashlib.sha1(seq).hexdigest()

    seqhash = pandas.Series(data=seqhashes, name='seqhash')
    seqhash.index.name = 'seqname'
    seq_info = seq_info.join(seqhash)

    group_by = ['seqhash']
    if args.group_by:
        group_by.extend(args.group_by.split(','))

    def choose_rep(df):
        if args.prefer_columns:
            df = df.sort_values(by=args.prefer_columns.split(','))
        rep = df.tail(1)
        rep['weight'] = len(df)
        return rep

    log.info('choosing seq_info representatives')
    seq_info = seq_info.groupby(
        by=group_by, group_keys=False).apply(choose_rep)
    seq_info = seq_info.drop('seqhash', axis=1)

    log.info('writing seqinfo')
    seq_info.to_csv(args.out_info, quoting=csv.QUOTE_NONNUMERIC)

    log.info('writing dedup file')
    with util.file_opener()(args.sequences) as sequences_in, \
            util.file_opener('w')(args.out) as sequences_out:
        for record in util.Counter(SeqIO.parse(sequences_in, 'fasta')):
            if record.name in seq_info.index:
                fasta_out = '>{}\n{}\n'.format(record.name, str(record.seq))
                sequences_out.write(fasta_out)
