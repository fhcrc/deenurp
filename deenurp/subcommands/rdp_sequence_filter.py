"""
Splits RDP sequence file into named/unnamed sections, filters by length, percent ambiguity.
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
    p.add_argument('fasta_file', help="""sequence file""", type=file_opener('r'))
    p.add_argument('seqinfo_file', help="""Sequence metadata""", type=file_opener('r'))
    p.add_argument('named_base')
    p.add_argument('unnamed_base')

    flt = p.add_argument_group('Filtering options')
    flt.add_argument('-a', '--prop-ambig-cutoff', default=0.01, type=float,
            help="""Maximum proportion of characters in sequence which may be
            ambiguous [default: %(default).2f]""")
    flt.add_argument('-l', '--min-length', type=int, help="""Minimum sequence
            length [default: %(default)d]""", default=1200)

def action(a):
    """
    Run
    """
    with a.fasta_file as fasta_fp, a.seqinfo_file as seqinfo_fp:
        sequences = Counter(SeqIO.parse(fasta_fp, 'fasta'))
        reader = csv.DictReader(seqinfo_fp)
        with open(a.named_base + '.fasta', 'w') as named_fa_fp, \
             open(a.named_base + '.seq_info.csv', 'w') as named_si_fp, \
             open(a.unnamed_base + '.fasta', 'w') as unnamed_fa_fp, \
             open(a.unnamed_base + '.seq_info.csv', 'w') as unnamed_si_fp:
            unnamed_writer = csv.DictWriter(unnamed_si_fp, reader.fieldnames,
                    quoting=csv.QUOTE_NONNUMERIC)
            unnamed_writer.writeheader()
            named_writer = csv.DictWriter(named_si_fp, reader.fieldnames,
                    quoting=csv.QUOTE_NONNUMERIC)
            named_writer.writeheader()

            accepted = 0
            rejected = 0
            for sequence, info in itertools.izip_longest(sequences, reader):
                # Check quality
                l = len(sequence)
                ambig_prop = count_ambiguous(str(sequence.seq)) / float(l)
                if ambig_prop > a.prop_ambig_cutoff or l < a.min_length:
                    rejected += 1
                    continue
                accepted += 1

                if info['tax_id'] != '':
                    # named
                    fa_fp = named_fa_fp
                    w = named_writer
                else:
                    fa_fp = unnamed_fa_fp
                    w = unnamed_writer
                w.writerow(info)
                SeqIO.write([sequence], fa_fp, 'fasta')

            logging.info("%d/%d passed (%.2f%%)", accepted, accepted+rejected,
                    accepted/float(accepted+rejected) * 100.0)
