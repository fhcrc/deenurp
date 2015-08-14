"""Splits RDP sequence file into named/unnamed sections,
filters by length, percent ambiguity.

Assumes sequences in fasta_file and records in seqinfo_file are in the
same order.

"""
import csv
import itertools
import logging

from Bio import SeqIO

from deenurp.util import Counter, file_opener


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def build_parser(p):
    p.add_argument(
        'fasta_file', help="""sequence file""", type=file_opener('r'))
    p.add_argument(
        'seqinfo_file', help="""Sequence metadata""", type=file_opener('r'))
    p.add_argument(
        '--named-seqs',
        default='named.seqs.fasta', help='[default %(default)s]')
    p.add_argument(
        '--named-info',
        default='named.seq_info.csv', help='[default %(default)s]')
    p.add_argument(
        '--unnamed-seqs',
        default='unnamed.seqs.fasta', help='[default %(default)s]')
    p.add_argument(
        '--unnamed-info',
        default='unnamed.seq_info.csv', help='[default %(default)s]')

    flt = p.add_argument_group('Filtering options')
    flt.add_argument('-a', '--prop-ambig-cutoff', default=0.01, type=float,
                     help="""Maximum proportion of characters in sequence which may be
                             ambiguous [default: %(default).2f]""")
    flt.add_argument('-l', '--min-length', type=int, help="""Minimum sequence
            length [default: %(default)d]""", default=1200)


def action(a):
    with a.fasta_file as fasta_fp, a.seqinfo_file as seqinfo_fp:
        sequences = Counter(SeqIO.parse(fasta_fp, 'fasta'))
        reader = csv.DictReader(seqinfo_fp)
        with open(a.named_seqs, 'w') as named_fa_fp, \
                open(a.named_info, 'w') as named_si_fp, \
                open(a.unnamed_seqs, 'w') as unnamed_fa_fp, \
                open(a.unnamed_info, 'w') as unnamed_si_fp:

            unnamed_writer = csv.DictWriter(
                unnamed_si_fp, reader.fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            unnamed_writer.writeheader()
            named_writer = csv.DictWriter(
                named_si_fp, reader.fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            named_writer.writeheader()

            accepted = 0
            rejected = 0
            for sequence, info in itertools.izip(sequences, reader):
                assert sequence.id == info['seqname']

                # Check quality
                l = len(sequence)
                ambig_prop = count_ambiguous(str(sequence.seq)) / float(l)
                if ambig_prop > a.prop_ambig_cutoff or l < a.min_length:
                    rejected += 1
                    continue
                accepted += 1

                if info['taxid_classified'] == 'True':
                    # named
                    fa_fp = named_fa_fp
                    w = named_writer
                else:
                    fa_fp = unnamed_fa_fp
                    w = unnamed_writer
                w.writerow(info)
                SeqIO.write([sequence], fa_fp, 'fasta')

            logging.info("%d/%d passed (%.2f%%)",
                         accepted,
                         accepted + rejected,
                         accepted / float(accepted + rejected) * 100.0)
