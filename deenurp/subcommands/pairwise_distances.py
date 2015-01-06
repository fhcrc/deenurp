"""Calculate a matrix of pairwise distances from unaligned sequences.

Output is a csv file with column labels corresponding to names in `seqs`
"""

from os.path import basename, splitext

import numpy

from deenurp.subcommands import filter_outliers
from deenurp.wrap import VSEARCH_IDDEF


def build_parser(parser):

    parser.add_argument('seqs', help="Fasta file")
    parser.add_argument('distmat', help="""File name for output
    distance matrix (csv with column labels)""")
    parser.add_argument(
        '-a', '--aligner', default='cmalign', choices=['cmalign', 'muscle', 'vsearch'])
    parser.add_argument(
        '--iddef', default=VSEARCH_IDDEF, type=int, choices=[0, 1, 2, 3, 4],
        help='vsearch: method for calculating pairwise identity [%(default)s]')
    parser.add_argument(
        '--threads', type=int, default=1,
        help="""number of threads (cmalign and vsearch only)""")


def action(args):

    pfx = splitext(basename(args.seqs))[0]

    if args.aligner == 'vsearch':
        taxa, distmat = filter_outliers.distmat_pairwise(
            args.seqs, pfx, args.aligner, iddef=args.iddef, threads=args.threads)
    elif args.aligner == 'cmalign':
        taxa, distmat = filter_outliers.distmat_cmalign(args.seqs, pfx, cpu=args.threads)
    elif args.aligner == 'muscle':
        taxa, distmat = filter_outliers.distmat_muscle(args.seqs, pfx)

    if args.distmat:
        numpy.savetxt(args.distmat, distmat, delimiter=',', header=','.join(taxa))
