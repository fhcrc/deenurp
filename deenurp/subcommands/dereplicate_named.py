"""Dereplicate reference sequences by clustering

"""

import logging
import os
import sys

import pandas as pd

from deenurp import uclust, wrap, util

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('seqs', help="""Named sequences""")
    parser.add_argument('seq_info', help="""Sequence info file""")
    parser.add_argument('-t', '--taxonomy',
                        help="""Taxonomy as taxtable; optional
                        if a grouping term is available in seq_info""")

    parser.add_argument('--seq-info-out',
                        help='subset of original seq_info')
    parser.add_argument('--derep-map-out',
                        help=('mapping of input sequences to dereplicated '
                              'representatives. `group` column corresponds to '
                              'the field identified by --group-on.'))
    parser.add_argument('--seqs-out',
                        default=sys.stdout,
                        type=util.file_opener('w'))

    parser.add_argument('-g', '--group-on', default='species',
                        help='Field in seq_info on which to group sequences')
    parser.add_argument('--id', default=1.0,
                        type=float, help="""Clustering identity between 0 and 1
                        [default: %(default).3f]""")
    parser.add_argument('-i', '--include', type=util.file_opener('r'),
                        help=('Optional file containing list '
                              'of group labels to include'))
    parser.add_argument('--threads',
                        help=('Number of threads to use for clustering each '
                              'group [default is one thread per '
                              'available CPU core]'))


def mocked_cluster_output(seqnames):
    columns = ['type', 'query_label', 'target_label']
    d = {'type': ['S' for s in seqnames],
         'query_label': seqnames, 'target_label': seqnames}
    return pd.DataFrame(d, columns=columns)


def cluster(seqfile, seqnames, identity=1.0, prefix='cluster-', threads=None):
    prefix = prefix.replace('/', '\\')  # / confuses the filesystem
    with util.ntf(prefix=prefix, suffix='.fasta') as fa, \
            util.ntf(prefix=prefix, suffix='.uc') as uc:
        wrap.esl_sfetch(seqfile, seqnames, fa)
        fa.flush()
        uclust.cluster(fa.name,
                       uc.name,
                       pct_id=identity,
                       pre_sorted=False,
                       quiet=True,
                       threads=threads)
        df = uclust.parse_uclust_as_df(uc)
        df = df[df.type != 'C']
        df = df[['type', 'query_label', 'target_label']]

        return df


def action(args):
    fa_idx = wrap.read_seq_file(args.seqs)

    dtype = {'gi': str, 'tax_id': str, 'species': str}
    seq_info = pd.read_csv(args.seq_info, dtype=dtype)
    info_cols = seq_info.columns

    if args.include:
        include = args.include.read().split()
        seq_info = seq_info.loc[seq_info[args.group_on].isin(include)]

    # join with taxonomy if provided
    if args.taxonomy:
        tax = pd.read_csv(args.taxonomy, dtype=str,).set_index('tax_id')
        seq_info = seq_info.join(tax, on='tax_id')

    grouped = seq_info.groupby(args.group_on, sort=False)

    frames = []
    for key, grp in grouped:
        # don't cluster groups represented by only one seq
        if grp.shape[0] == 1:
            clusters = mocked_cluster_output(grp['seqname'])
        else:
            # TODO: is longer necessarily better?
            by = []
            ascending = []
            for c, o in [['is_type', False],
                         ['ambig_count', True],
                         ['length', False]]:
                if c in grp.columns:
                    by.append(c)
                    ascending.append(o)
            if by:
                grp = grp.sort_values(
                    by=by,
                    ascending=ascending)
                clusters = cluster(
                    args.seqs,
                    grp['seqname'],
                    fa_idx,
                    identity=args.id,
                    prefix='{}-'.format(key),
                    threads=args.threads)

        clusters['group'] = key
        frames.append(clusters)

    all_clusters = pd.concat(frames)

    if args.derep_map_out:
        all_clusters.columns = ['type', 'seqname', 'seed', 'group']
        all_clusters.to_csv(args.derep_map_out, header=True, index=False)

    if args.seq_info_out:
        seq_info = seq_info[seq_info['seqname'].isin(all_clusters['seed'])]
        seq_info.to_csv(args.seq_info_out, columns=info_cols, index=False)

    wrap.esl_sfetch(
        args.seqs,
        all_clusters['seed'].unique(),
        args.seqs_out,
        fa_idx)
