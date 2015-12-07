"""Align nucleotide sequences using cmalign

"""

import logging

from Bio import SeqIO

from .. import wrap, util

log = logging.getLogger(__name__)


def cmalign(infile, outfile, cpu):

    with util.ntf(suffix='.sto') as a_sto, open(outfile, 'w') as a_fasta:
        scores = wrap.cmalign_files(infile, a_sto.name, cpu=cpu)
        SeqIO.convert(a_sto, 'stockholm', a_fasta, 'fasta')
        a_fasta.flush()

    return scores


def build_parser(p):
    p.add_argument('infile', help='input file in fasta format')
    p.add_argument('outfile', help='output file in fasta format')
    p.add_argument('--scores',
                   help='csv file containing cmalign scores')
    p.add_argument('-t', '--threads', type=int, default=wrap.CMALIGN_THREADS,
                   help='number of processes [%(default)s]')


def action(a):
    scores = cmalign(a.infile, a.outfile, cpu=a.threads)
    if a.scores:
        scores.to_csv(a.scores)
