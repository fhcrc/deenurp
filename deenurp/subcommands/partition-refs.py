"""Splits sequence files into partions and optionally
filters by length and percent ambiguity.

partions: named, unnamed, types, published

"""
from Bio import SeqIO

from deenurp.util import Counter, file_opener, read_csv


def build_parser(p):
    # inputs
    p.add_argument(
        'fasta',
        metavar='FASTA',
        help="""sequence file""",
        type=file_opener('r'))
    p.add_argument(
        'seqinfo',
        metavar='CSV',
        help="""Sequence metadata""")
    p.add_argument(
        '--references',
        metavar='CSV',
        help='csv file with columns: version,pubmed_id')

    # outputs
    info_outs = p.add_argument_group('outputs for seq_info\'s')
    info_outs.add_argument(
        '--named-info',
        metavar='CSV',
        type=file_opener('w'),
        help='taxid_classified column is True')
    info_outs.add_argument(
        '--unnamed-info',
        metavar='CSV',
        type=file_opener('w'),
        help='taxid_classified column is False')
    info_outs.add_argument(
        '--type-info',
        type=file_opener('w'),
        metavar='CSV',
        help='rows where is_type column is True')
    info_outs.add_argument(
        '--published-info',
        type=file_opener('w'),
        metavar='CSV',
        help="""requires references.csv. Return seq_info with pubmed_ids""")
    info_outs.add_argument(
        '--references-info',
        type=file_opener('w'),
        metavar='CSV',
        help=('requires references.csv. '
              'Return columns [version, accession, pubmed_id]'))

    seq_outs = p.add_argument_group('outputs for sequences')
    seq_outs.add_argument(
        '--named-seqs',
        type=file_opener('w'),
        metavar='FASTA',
        help='where taxid_classified column is True')
    seq_outs.add_argument(
        '--unnamed-seqs',
        type=file_opener('w'),
        metavar='FASTA',
        help='where taxid_classified column is False')
    seq_outs.add_argument(
        '--type-seqs',
        type=file_opener('w'),
        metavar='FASTA',
        help='where is_type column is True')
    seq_outs.add_argument(
        '--published-seqs',
        type=file_opener('w'),
        metavar='FASTA',
        help="""requires references.csv. Return sequences with pubmed_ids""")

    # filtering switches
    flt = p.add_argument_group('filtering options')
    flt.add_argument('-a', '--prop-ambig-cutoff',
                     type=float,
                     help=('Maximum proportion of characters in '
                           'sequence which may be ambiguous'))
    flt.add_argument('-l', '--min-length',
                     type=int,
                     help='Minimum sequence length')


def action(args):
    # load seq_info
    seq_info = read_csv(args.seqinfo)
    seq_info = seq_info.set_index('id')

    # raw min_length filtering
    if args.min_length:
        seq_info = seq_info[seq_info['length'] > args.min_length]

    # raw prop_ambig filtering
    if args.prop_ambig_cutoff:
        seq_info['prop_ambig'] = seq_info['ambig_count'] / seq_info['length']
        seq_info = seq_info[seq_info['prop_ambig'] < args.prop_ambig_cutoff]
        seq_info = seq_info.drop('prop_ambig', axis=1)

    if args.unnamed_seqs or args.unnamed_info:
        unnamed_info = seq_info[~seq_info['taxid_classified']]
        if args.unnamed_info:
            unnamed_info.to_csv(args.unnamed_info)

    if args.named_seqs or args.named_info:
        named_info = seq_info[seq_info['taxid_classified']]
        if args.named_info:
            named_info.to_csv(args.named_info)

    if args.type_seqs or args.type_info:
        type_info = seq_info[seq_info['is_type']]
        if args.type_info:
            type_info.to_csv(args.type_info)

    if (args.published_seqs or args.references_info) and args.references:
        references = read_csv(args.references,
                              usecols=['version', 'pubmed_id'])
        published_info = seq_info.reset_index()
        published_info = published_info.merge(
            references, on='version', how='left')
        published_info = published_info[published_info['pubmed_id'].notnull()]

        if args.references_info:
            ref_cols = references.columns.tolist() + ['accession']
            references_info = published_info[ref_cols].drop_duplicates()
            references_info.to_csv(args.references_info, index=False)

        published_info = published_info.drop('pubmed_id', axis=1)
        published_info = published_info.drop_duplicates().set_index('id')
        if args.published_info:
            published_info.to_csv(args.published_info)

    if (args.unnamed_seqs or args.named_seqs or
            args.type_seqs or (args.published_seqs and args.references)):
        for fa in Counter(SeqIO.parse(args.fasta, 'fasta')):
            if args.unnamed_seqs and fa.id in unnamed_info.index:
                SeqIO.write([fa], args.unnamed_seqs, 'fasta')

            if args.named_seqs and fa.id in named_info.index:
                SeqIO.write([fa], args.named_seqs, 'fasta')

            if args.type_seqs and fa.id in type_info.index:
                SeqIO.write([fa], args.type_seqs, 'fasta')

            if (args.published_seqs and
                    args.references and
                    fa.id in published_info.index):
                SeqIO.write([fa], args.published_seqs, 'fasta')
