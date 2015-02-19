"""Remove mis-annotated or malformed sequences from a reference database.

Inputs are a database of reference sequences in fasta format, sequence
annotation, and a taxonomy file.

Sequences are first grouped at a specified taxonomic rank (the default
is species). For each group, all pairwise distances are calculated
using a choice of alignment strategies. If ``--aligner=cmalign``, a
multiple alignment is created using a 16S rRNA alignment profile. A
multiple alignment of non-16S sequences can be created using
``--aligner=muscle``, but note that alignment of large numbers of
sequences may be slow. For both multiple alignment strategies,
pairwise distances are calculated using ``FastTree
-makematrix``. Alternatively, ``--aligner=vsearch`` will calculate
pairwise distances using global pairwise alignments.

Given a pairwise distance matrix, two strategies are available for
outlier detection:

* ``--strategy=radius`` finds the medoid of the group and discards all
  other sequences with a distance from the medoid greater than some
  threshold.
* ``--strategy=cluster`` performs single-linkage hierarchical
  clustering at a specified threshold T. All clusters of size 1 are
  discarded. In addition, medoids of each remaining cluster are
  identified, and all members of any cluster whose medoid has a
  distance of 1.5 * T from the medoid of the largest cluster is also
  discarded.

There are two options for defining the threshold T:

* By default, T is defined as a percentile of the distribution of
  distances from the group medoid to all the other sequences in the
  group (``--distance-percentile``). Bounds on this value may be
  defined using ``--{min,max}-distance``.
* Alternatively, a fixed value may be defined using ``--distance-cutoff``

At the species level, T will typically be 0.01 to 0.02; 0.015 is a
reasonable to use as a fixed value.

A filtered subset of the sequences are provided as output. The file
specified by ``--filtered-seqinfo`` contains the filtered subset of
``seqinfo_file``. ``--detailed-seqinfo`` contains the following additional fields:
* centroid - the group centroid
* dist - distance to the group centroid
* is_out - whether the sequence is an outlier
* {rank}_id - the tax_id at the taxonomic rank used for grouping
* x, y - coordinates calculated by multidimensional scaling of the
  pairwise distance matrix.

"""

import argparse
import os
import os.path
import sys
import shutil
import logging

log = logging

from concurrent import futures
import numpy
import pandas as pd

from Bio import SeqIO
from taxtastic.taxtable import TaxNode

from deenurp import config, wrap, util, outliers
from deenurp.wrap import (CMALIGN_THREADS, MUSCLE_MAXITERS,
                          VSEARCH, VSEARCH_IDDEF, VSEARCH_THREADS)

DEFAULT_RANK = 'species'
DEFAULT_ALIGNER = 'cmalign'
DROP = 'drop'
KEEP = 'keep'
BLAST6NAMES = ['query', 'target', 'pct_id', 'align_len', 'mismatches', 'gaps',
               'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bit_score']


def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""")
    p.add_argument(
        'taxonomy', help="""Taxtable""", type=argparse.FileType('r'))

    output_group = p.add_argument_group('output options')
    output_group.add_argument(
        'output_fp', help="""Destination for sequences""",
        type=argparse.FileType('w'))
    output_group.add_argument(
        '--filtered-seqinfo', type=argparse.FileType('w'),
        help="""Path to write filtered sequence info""")
    output_group.add_argument(
        '--detailed-seqinfo', type=argparse.FileType('w'),
        help="""Sequence info, including filtering details""")

    filter_group = p.add_argument_group('filtering options')
    filter_group.add_argument('--filter-rank', default=DEFAULT_RANK)
    filter_group.add_argument(
        '--strategy', default='radius', choices=['radius', 'cluster'],
        help="""Strategy for outlier detection.""")
    filter_group.add_argument(
        '--distance-percentile', type=float, default=90.0,
        help="""Define distance cutoff as a percentile of the
        distribution of distances from medoid to others [default:
        %(default)s]""")
    filter_group.add_argument(
        '--min-distance', type=float, default=0.01,
        help="""Minimum distance when calculating as a percentile
        [default: %(default)s]""")
    filter_group.add_argument(
        '--max-distance', type=float, default=0.10,
        help="""Maximum distance when calculating as a percentile
        [default: %(default)s]""")
    filter_group.add_argument(
        '--distance-cutoff', type=float,
        help="""Distance threshold from cluster centroid
        (--strategy='radius') or distance parameter for hierarchical
        clustering (--strategy='cluster'). Overrides distance
        calculation using --distance-percentile if provided [default:
        %(default)s]""")

    aligner_group = p.add_argument_group("aligner-specific options")
    aligner_group.add_argument(
        '--aligner', help='multiple alignment tool [%(default)s]',
        default=DEFAULT_ALIGNER, choices=['cmalign', 'muscle', 'vsearch'])
    aligner_group.add_argument(
        '--executable',
        help='Optional absolute or relative path to the alignment tool executable')
    # aligner_group.add_argument(
    #     '--maxiters', default=MUSCLE_MAXITERS, type=int,
    #     help='muscle: value for -maxiters [%(default)s]')
    # aligner_group.add_argument(
    #     '--iddef', default=VSEARCH_IDDEF, type=int, choices=[0, 1, 2, 3, 4],
    #     help='vsearch: method for calculating pairwise identity [%(default)s]')

    rare_group = p.add_argument_group("rare taxa")
    rare_group.add_argument(
        '--min-seqs-for-filtering', type=int, default=5,
        help="""Minimum number of sequences perform distance-based
            medoid-filtering on [default: %(default)d]""")
    rare_group.add_argument(
        '--rare-taxon-action', choices=(KEEP, DROP), default=KEEP,
        help="""Action to perform when a taxon has <
            '--min-seqs-to-filter' representatives. [default: %(default)s]""")

    p.add_argument('--threads', type=int, default=config.DEFAULT_THREADS,
                   help="""number of taxa to process concurrently (one
                   process per multiple alignment)""")


def sequences_above_rank(taxonomy, rank=DEFAULT_RANK):
    """
    Generate the sequence ids from taxonomy whose rank is above specified rank.
    """
    ranks = taxonomy.ranks
    r_index = ranks.index(rank)
    assert r_index >= 0

    def above_rank(node):
        n_index = ranks.index(node.rank)
        assert n_index >= 0
        return n_index < r_index

    for n in taxonomy:
        if above_rank(n):
            for sequence_id in n.sequence_ids:
                yield sequence_id


def distmat_muscle(sequence_file, prefix, maxiters=MUSCLE_MAXITERS):

    with util.ntf(prefix=prefix, suffix='.fasta') as a_fasta:
        wrap.muscle_files(sequence_file, a_fasta.name, maxiters=maxiters)
        a_fasta.flush()

        taxa, distmat = outliers.fasttree_dists(a_fasta.name)

    return taxa, distmat


def distmat_cmalign(sequence_file, prefix, cpu=CMALIGN_THREADS):

    with util.ntf(prefix=prefix, suffix='.aln') as a_sto, \
            util.ntf(prefix=prefix, suffix='.fasta') as a_fasta:

        wrap.cmalign_files(sequence_file, a_sto.name, cpu=cpu)
        # FastTree requires FASTA
        SeqIO.convert(a_sto, 'stockholm', a_fasta, 'fasta')
        a_fasta.flush()

        taxa, distmat = outliers.fasttree_dists(a_fasta.name)

    return taxa, distmat


class UsearchError(Exception):
    pass


def parse_usearch_allpairs(filename, seqnames):
    """Read output of ``usearch -allpairs_global -blast6out`` and return a
    square distance matrix. ``seqnames`` determines the marginal order
    of sequences in the matrix.

    """

    data = pd.io.parsers.read_table(filename, header=None, names=BLAST6NAMES)
    data['dist'] = pd.Series(1.0 - data['pct_id'] / 100.0, index=data.index)

    # for each sequence pair, select the longest alignment if there is
    # more than one (chooses first occurrence if there are two the same length).
    maxidx = data.groupby(['query', 'target']).apply(lambda x: x['align_len'].idxmax())
    data = data.iloc[maxidx]

    if set(seqnames) != set(data['query']) | set(data['target']):
        # shutil.copy(filename, '.')
        raise UsearchError(
            'some sequences are missing from the output ({})'.format(filename))

    nseqs = len(seqnames)
    distmat = numpy.repeat(0.0, nseqs ** 2)
    distmat.shape = (nseqs, nseqs)
    ii = pd.match(data['query'], seqnames)
    jj = pd.match(data['target'], seqnames)

    # usearch_allpairs_files returns comparisons corresponding to a
    # triangular matrix, whereas vsearch_allpairs_files returns all
    # comparisons. Here we convert both to a square matrix.
    if data.shape[0] == nseqs * nseqs:
        distmat[ii, jj] = data['dist']
    elif data.shape[0] == (nseqs * (nseqs - 1)) / 2:
        distmat[ii, jj] = data['dist']
        distmat[jj, ii] = data['dist']
    else:
        raise UsearchError(
            'not all pairwise comparisons are represented ({})'.format(filename))

    return distmat


def distmat_pairwise(sequence_file, prefix, aligner,
                     executable=None, iddef=VSEARCH_IDDEF, threads=VSEARCH_THREADS):
    """Calculate a distance matrix using either usearch or vsearch.

    """

    assert aligner in {'vsearch'}, 'invalid aligner'
    executable = executable or {'vsearch': VSEARCH}[aligner]

    with open(sequence_file) as sf, util.ntf(
            prefix=prefix, suffix='.blast6out') as uc:

        if aligner == 'vsearch':
            wrap.vsearch_allpairs_files(sequence_file, uc.name, executable,
                                        iddef=iddef, threads=threads)

        uc.flush()

        taxa = [seq.id for seq in SeqIO.parse(sf, 'fasta')]
        try:
            distmat = parse_usearch_allpairs(uc.name, taxa)
        except UsearchError, e:
            logging.error(e)
            shutil.copy(sequence_file, '.')
            raise UsearchError

    return taxa, distmat


def filter_sequences(tax_id,
                     sequence_file=None,
                     distmat=None,
                     taxa=None,
                     strategy='radius',
                     cutoff=None,
                     percentile=None,
                     min_radius=0.0,
                     max_radius=None,
                     cluster_type='single',
                     aligner='cmalign',
                     executable=None,
                     maxiters=MUSCLE_MAXITERS, iddef=VSEARCH_IDDEF):
    """
    Return a list of sequence names identifying outliers.
    """

    assert aligner in {'cmalign', 'muscle', 'vsearch'}, 'invalid aligner'
    assert strategy in {'radius', 'cluster'}, 'invalid strategy'

    prefix = '{}_'.format(tax_id)

    if distmat is None:
        if aligner == 'cmalign':
            taxa, distmat = distmat_cmalign(sequence_file, prefix)
        elif aligner == 'muscle':
            taxa, distmat = distmat_muscle(sequence_file, prefix, maxiters)
        elif aligner == 'vsearch':
            taxa, distmat = distmat_pairwise(
                sequence_file, prefix, aligner, executable, iddef)
    else:
        assert taxa is not None

    if cutoff is not None:
        cutoff = cutoff
    elif percentile is not None:
        cutoff = outliers.scaled_radius(distmat, percentile, min_radius, max_radius)
    else:
        raise ValueError('must provide either cutoff or percentile')

    log.info('cutoff={} ({})'.format(
        cutoff, 'calculated' if percentile else 'pre-defined'))

    log.info('strategy: {}'.format(strategy))
    if strategy == 'radius':
        medoid, dists, is_out = outliers.outliers(distmat, cutoff)
    elif strategy == 'cluster':
        medoid, dists, is_out = outliers.outliers_by_cluster(
            distmat, t=cutoff, D=cutoff * 1.5,
            min_size=2, cluster_type=cluster_type)

    assert len(is_out) == len(taxa)

    result = pd.DataFrame({
        'seqname': taxa,
        'centroid': numpy.repeat(taxa[medoid], len(taxa)),
        'dist': dists,
        'is_out': is_out})

    mds = outliers.mds(distmat, taxa)
    result = pd.merge(result, mds, how='left', on='seqname')

    return result


def mock_filter(seqs, keep):
    """Return a DataFrame with the same structure as the output of
    `filter_worker`. The value of the 'is_out' column is (not keep)
    for all seqs.

    """

    empty = numpy.repeat(numpy.nan, len(seqs))
    return pd.DataFrame({
        'seqname': seqs,
        'centroid': empty,
        'dist': empty,
        'is_out': numpy.repeat(not keep, len(seqs))})


def filter_worker(tax_id,
                  sequence_file,
                  seqs,
                  strategy,
                  distance_cutoff,
                  percentile,
                  min_radius,
                  max_radius,
                  cluster_type,
                  aligner,
                  executable):
    """
    Worker task for running filtering tasks.

    Arguments:
    :tax_id: string identifying tax_id of these sequences
    :sequence_file: Complete sequence file
    :seqs: sequence names representing this tax_id
    :strategy: 'radius' or 'cluster'
    :distance_cutoff: Distance cutoff for medoid filtering
    :percentile: Used to calculate distance cutoff (``distance_cutoff`` cannot also be provided)
    :{min, max}_radius: Bounds of calculated distance cutoff
    :cluster_type: clustering algorithm (method of scipy.cluster.hierarchy)
    :aligner: name of alignment program
    :executable: name or path of executable for alignmment program

    :returns: output of ``filter_sequences()``
    """

    prefix = '{}_'.format(tax_id)

    with util.ntf(prefix=prefix, suffix='.fasta') as tf:
        # Extract sequences
        wrap.esl_sfetch(sequence_file, seqs, tf)
        tf.flush()

        filtered = filter_sequences(
            tax_id,
            sequence_file=tf.name,
            strategy=strategy,
            cutoff=distance_cutoff,
            percentile=percentile,
            min_radius=min_radius,
            max_radius=max_radius,
            aligner=aligner,
            executable=executable)

        return filtered


def action(a):
    # remove .ssi index for sequence file if it exists
    try:
        os.remove(a.sequence_file + '.ssi')
    except OSError:
        pass

    # Load taxonomy
    with a.taxonomy as fp:
        taxonomy = TaxNode.from_taxtable(fp)
        logging.info('Loaded taxonomy')

    # Load sequences into taxonomy
    with open(a.seqinfo_file) as fp:
        taxonomy.populate_from_seqinfo(fp)
        logging.info('Added %d sequences',
                     sum(1 for i in taxonomy.subtree_sequence_ids()))

    executable = a.executable or {'cmalign': 'cmalign',
                                  'muscle': 'muscle',
                                  'vsearch': wrap.VSEARCH}[a.aligner]

    outcomes = []  # accumulate DatFrame objects

    filter_rank_col = '{}_id'.format(a.filter_rank)

    # Sequences which are classified above the desired rank should just be kept
    names_above_rank = list(sequences_above_rank(taxonomy, a.filter_rank))
    logging.info('Keeping %d sequences classified above %s',
                 len(names_above_rank), a.filter_rank)
    above_rank = mock_filter(names_above_rank, keep=True)
    above_rank[filter_rank_col] = pd.Series(numpy.nan, index=above_rank.index)
    outcomes.append(above_rank)

    # For each filter-rank, filter
    nodes = [i for i in taxonomy if i.rank == a.filter_rank]

    # Filter each tax_id, running ``--threads`` tasks in parallel
    with futures.ThreadPoolExecutor(a.threads) as executor:
        # dispatch a pool of tasks
        futs = {}
        for i, node in enumerate(nodes):
            seqs = frozenset(node.subtree_sequence_ids())

            if not seqs:
                logging.debug("No sequences for %s (%s)", node.tax_id, node.name)
                continue
            elif len(seqs) < a.min_seqs_for_filtering:
                logging.debug('%d sequence(s) for %s (%s) [action: %s]',
                              len(seqs), node.tax_id, node.name, a.rare_taxon_action)
                f = executor.submit(mock_filter, seqs=list(seqs),
                                    keep=a.rare_taxon_action == KEEP)
            else:
                # `f` is a DataFrame (output of filter_sequences)
                f = executor.submit(
                    filter_worker,
                    tax_id=node.tax_id,
                    sequence_file=a.sequence_file,
                    seqs=seqs,
                    strategy=a.strategy,
                    distance_cutoff=a.distance_cutoff,
                    percentile=a.distance_percentile,
                    min_radius=a.min_distance,
                    max_radius=a.max_distance,
                    cluster_type='single',
                    aligner=a.aligner,
                    executable=executable)

            futs[f] = {'n_seqs': len(seqs), 'node': node}

        # log results for each tax_id as tasks complete
        complete = 0
        while futs:
            done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
            complete += len(done)
            sys.stderr.write('{0:8d}/{1:8d} taxa completed\r'.format(
                complete, complete + len(pending)))
            for f in done:
                if f.exception():
                    logging.exception("Error in child process: %s", f.exception())
                    executor.shutdown(False)
                    raise f.exception()

                info = futs.pop(f)
                filtered = f.result()  # here's the DataFrame again...

                # add a column for tax_d at filter_rank
                filtered[filter_rank_col] = pd.Series(
                    info['node'].tax_id, index=filtered.index)
                outcomes.append(filtered)

                kept = frozenset(filtered.seqname[~filtered.is_out])
                if len(kept) == 0:
                    logging.info('Pruned all %d sequences for %s (%s)',
                                 info['n_seqs'], info['node'].tax_id,
                                 info['node'].name)
                elif len(kept) != info['n_seqs']:
                    logging.info('Pruned %d/%d sequences for %s (%s)',
                                 info['n_seqs'] - len(kept), info['n_seqs'],
                                 info['node'].tax_id, info['node'].name)

    all_outcomes = pd.concat(outcomes, ignore_index=True)
    all_outcomes.set_index('seqname', inplace=True)

    # all input sequences should be in the output
    assert {s for node in taxonomy for s in node.sequence_ids} == set(all_outcomes.index)

    kept_ids = set(all_outcomes.index[~all_outcomes.is_out])

    with a.output_fp as fp:
        # Extract all of the sequences that passed.
        logging.info('Extracting %d sequences', len(kept_ids))
        wrap.esl_sfetch(a.sequence_file, kept_ids, fp)

    # Filter seqinfo for sequences that passed.
    seqinfo = pd.io.parsers.read_csv(
        a.seqinfo_file, dtype={'seqname': str, 'tax_id': str})
    seqinfo.set_index('seqname', inplace=True)

    merged = seqinfo.join(all_outcomes, lsuffix='.left')

    # csv output
    if a.filtered_seqinfo:
        merged[~merged.is_out].to_csv(a.filtered_seqinfo, columns=seqinfo.columns)

    if a.detailed_seqinfo:
        merged.to_csv(a.detailed_seqinfo)
        # merged.to_csv(a.detailed_seqinfo[['seqname', 'is_out', 'centroid', 'dist', 'x', 'y']])
