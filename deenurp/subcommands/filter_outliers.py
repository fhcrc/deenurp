"""Remove mis-annotated or malformed sequences from a reference database.

algorithm
=========

Inputs are a database of reference sequences in fasta format, sequence
annotation, and a taxonomy file.

Sequences are first grouped at a specified taxonomic rank (the default
is species). For each group of sequences, all pairwise distances are
calculated using a choice of strategies. If ``--aligner=cmalign``, a
multiple alignment is created using a 16S rRNA alignment profile. A
multiple alignment of non-16S sequences can be created using
``--aligner=muscle``, but note that alignment of large numbers of
sequences may be slow. For both multiple alignment strategies,
pairwise distances are calculated using ``FastTree
-makematrix``. Alternatively, ``--aligner=vsearch`` will calculate
pairwise distances using global pairwise alignments. This tends to be
faster for small groups of sequences, and seems acceptable for up to
1000 or so records per group.

Given a pairwise distance matrix, two strategies are available for
outlier detection:

* ``--strategy=radius`` finds the medoid of the group and discards all
  other sequences with a distance from the medoid greater than some
  threshold.
* ``--strategy=cluster`` performs single-linkage hierarchical
  clustering at a specified threshold T. All clusters of size 1 are
  discarded. In addition, medoids of each remaining cluster are
  identified, and all members of any cluster whose medoid has a
  distance > 1.5 * T from the medoid of the largest cluster is also
  discarded.

There are two options for defining the threshold T:

* By default, T is defined as a percentile of the distribution of
  distances from the group medoid to all the other sequences in the
  group (``--distance-percentile``). Bounds on this value may be
  defined using ``--{min,max}-distance``.
* Alternatively, a fixed value may be defined using ``--distance-cutoff``

For full-length 16S rRNA sequences grouped at the species level, T
will typically be 0.01 to 0.02; 0.015 is a reasonable to use as a
fixed value.

A filtered subset of the sequences are provided as output. The file
specified by ``--filtered-seqinfo`` contains the filtered subset of
``seqinfo_file``. ``--detailed-seqinfo`` contains
the following additional fields:
* centroid - the group centroid
* dist - distance to the group centroid
* is_out - whether the sequence is an outlier
* {rank} - the tax_id at the taxonomic rank used for grouping
* x, y - coordinates calculated by multidimensional scaling of the
  pairwise distance matrix.
* cluster - an integer identifying a cluster of reads when
  ``--strategy=cluster`` is used, or the group centroid using
  ``--strategy=radius``

This can all take a while if there are many sequences. If a pool of
candidate sequences needs to be updated, use ``--previous-details`` to
provide the output of ``--detailed-seqinfo`` from a previous run to
avoid re-analyzing tax_ids represented by the same set of sequences.

runtime parameters
==================

Parallel execution is supported at two levels: ``--jobs`` defines the
number of sequence groups that are processed in parallel, and
``--threads-per-job`` determines the number of cpus or threads
allocated to each job (eg via ``vsearch --threads`` or ``cmalign
--cpu``). No effort is made to avoid exceeding available resources, so
the user should consider the product of these two parameters.

"""

import argparse
import numpy
import os
import pandas as pd
import sys
import shutil
import logging
import csv
import traceback

from Bio import SeqIO
from concurrent import futures

from taxtastic.taxtable import TaxNode as _TaxNode
from .. import config, wrap, util, outliers

log = logging.getLogger(__name__)

DEFAULT_RANK = 'species'
DEFAULT_ALIGNER = 'cmalign'
DROP = 'drop'
KEEP = 'keep'
BLAST6NAMES = ['query', 'target', 'pct_id', 'align_len', 'mismatches', 'gaps',
               'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bit_score']


# monkey patch class from taxtastic to warn when tax_id is missing
# from the taxonomy
# TODO: remove this and fix in taxtastic
class TaxNode(_TaxNode):

    def populate_from_seqinfo(self, seqinfo, seqnames):
        """Populate sequence_ids below this node from a seqinfo file object."""
        for row in csv.DictReader(seqinfo):
            if row['seqname'] in seqnames:
                node = self.index.get(row['tax_id'])
                if node:
                    node.sequence_ids.add(row['seqname'])
                else:
                    log.warning('tax_id {tax_id} (sequence {seqname}) '
                                'was not found in the taxonomy'.format(**row))


def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""")
    p.add_argument(
        'taxonomy', help="""Taxtable""", type=argparse.FileType('r'))

    input_group = p.add_argument_group('other inputs')
    input_group.add_argument(
        '--previous-details', metavar='FILE',
        help="""Output of --detailed-seqinfo from a previous run. If provided, use
        the previous results for any tax_ids represented by the same set of
        sequences in both sequence_file and this file.""")

    output_group = p.add_argument_group('output options')
    output_group.add_argument(
        '--output-seqs', help="""REQUIRED destination for sequences""",
        required=True,
        type=argparse.FileType('w'), metavar='FILE')
    output_group.add_argument(
        '--filtered-seqinfo', type=argparse.FileType('w'), metavar='FILE',
        help="""Path to write filtered sequence info""")
    output_group.add_argument(
        '--detailed-seqinfo', metavar='FILE',
        help="""Sequence info, including filtering details""")

    filter_group = p.add_argument_group('filtering options')

    # TODO: write more info in help about what this column is doing
    filter_group.add_argument(
        '--filter-rank', default=DEFAULT_RANK, help='[%(default)s]')
    filter_group.add_argument(
        '--no-filter',
        help=('file or list of sequence taxids '
              'to pass through without filtering'))
    filter_group.add_argument(
        '--strategy', default='radius', choices=['radius', 'cluster'],
        help="""Strategy for outlier detection. """)
    filter_group.add_argument(
        '--cluster-type', default='single',
        choices=['single', 'RobustSingleLinkage'],
        help="""Identifies the clustering algorithm for
        --strategy=cluster. 'single' for
        'scipy.cluster.hierarchy.single' or 'RobustSingleLinkage' for
        'hdbscan.RobustSingleLinkage'""")
    filter_group.add_argument(
        '--distance-percentile', type=float, default=50.0,
        help="""Define distance cutoff as a percentile of the
        distribution of distances from medoid to others [default:
        %(default)s]""")
    filter_group.add_argument(
        '--min-distance', type=float, default=0.01,
        help="""Minimum distance when calculating as a percentile
        [default: %(default)s]""")
    filter_group.add_argument(
        '--max-distance', type=float, default=0.05,
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
        help=('Optional absolute or relative path '
              'to the alignment tool executable'))

    rare_group = p.add_argument_group("rare taxa")
    rare_group.add_argument(
        '--min-seqs-for-filtering', type=int, default=5,
        help="""Minimum number of sequences perform distance-based
            medoid-filtering on [default: %(default)d]""")
    rare_group.add_argument(
        '--rare-taxon-action', choices=(KEEP, DROP), default=KEEP,
        help="""Action to perform when a taxon has <
            '--min-seqs-to-filter' representatives. [default: %(default)s]""")

    p.add_argument('-j', '--jobs', type=int, default=config.DEFAULT_THREADS,
                   help="""number of taxa to process concurrently [default %(default)s]""")
    p.add_argument('-t', '--threads-per-job', type=int, default=4,
                   help="""number of threads per job (eg, value to pass 'cmalign --cpu')
                   [default %(default)s]""")


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


def distmat_muscle(sequence_file, prefix, maxiters=wrap.MUSCLE_MAXITERS):

    with util.ntf(prefix=prefix, suffix='.fasta') as a_fasta:
        wrap.muscle_files(sequence_file, a_fasta.name, maxiters=maxiters)
        a_fasta.flush()

        taxa, distmat = outliers.fasttree_dists(a_fasta.name)

    return taxa, distmat


def distmat_cmalign(
        sequence_file,
        prefix,
        cpu=wrap.CMALIGN_THREADS,
        min_bitscore=10):

    with util.ntf(prefix=prefix, suffix='.aln') as a_sto, \
            util.ntf(prefix=prefix, suffix='.fasta') as a_fasta:

        scores = wrap.cmalign_files(sequence_file, a_sto.name, cpu=cpu)

        low_scores = scores['bit_sc'] < min_bitscore
        if low_scores.any():
            msg = 'The following sequences aligned with bit score < {}: {}'
            log.warning(msg.format(min_bitscore, scores[low_scores].index))

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

    data = pd.read_table(filename, header=None, names=BLAST6NAMES)
    data['dist'] = pd.Series(
        1.0 - data['pct_id'] / 100.0, index=data.index)

    # for each sequence pair, select the longest alignment if there is
    # more than one (chooses first occurrence if there are two the same
    # length).
    maxidx = data.groupby(['query', 'target']).apply(
        lambda x: x['align_len'].idxmax())
    data = data.iloc[maxidx]

    if set(seqnames) != set(data['query']) | set(data['target']):
        # shutil.copy(filename, '.')
        raise UsearchError(
            'some sequences are missing from the output ({})'.format(filename))

    nseqs = len(seqnames)
    distmat = numpy.repeat(0.0, nseqs ** 2)
    distmat.shape = (nseqs, nseqs)

    idx = dict(zip(seqnames, range(nseqs)))
    ii = [idx[name] for name in data['query']]
    jj = [idx[name] for name in data['target']]

    # usearch_allpairs_files returns comparisons corresponding to a
    # triangular matrix, whereas vsearch_allpairs_files returns all
    # comparisons. Here we convert both to a square matrix.
    if data.shape[0] == nseqs * nseqs:
        distmat[ii, jj] = data['dist']
    elif data.shape[0] == (nseqs * (nseqs - 1)) / 2:
        distmat[ii, jj] = data['dist']
        distmat[jj, ii] = data['dist']
    else:
        msg = 'not all pairwise comparisons are represented ({})'
        raise UsearchError(msg.format(filename))

    return distmat


def distmat_pairwise(sequence_file,
                     prefix, aligner,
                     executable=None,
                     iddef=wrap.VSEARCH_IDDEF,
                     threads=wrap.VSEARCH_THREADS):
    """Calculate a distance matrix using either usearch or vsearch.

    """

    assert aligner in {'vsearch'}, 'invalid aligner'
    executable = executable or {'vsearch': wrap.VSEARCH}[aligner]

    with open(sequence_file) as sf, util.ntf(
            prefix=prefix, suffix='.blast6out') as uc:

        if aligner == 'vsearch':
            wrap.vsearch_allpairs_files(sequence_file, uc.name, executable,
                                        iddef=iddef, threads=threads)

        uc.flush()

        taxa = [seq.id for seq in SeqIO.parse(sf, 'fasta')]
        try:
            distmat = parse_usearch_allpairs(uc.name, taxa)
        except UsearchError as e:
            log.error(e)
            shutil.copy(sequence_file, '.')
            raise UsearchError

    return taxa, distmat


def filter_sequences(tax_id,
                     sequence_file=None,
                     distmat=None,
                     taxa=None,
                     strategy='radius',
                     cluster_type='single',
                     cutoff=None,
                     percentile=None,
                     min_radius=0.0,
                     max_radius=None,
                     aligner='cmalign',
                     executable=None,
                     maxiters=wrap.MUSCLE_MAXITERS,
                     iddef=wrap.VSEARCH_IDDEF,
                     threads=None):
    """
    Return a list of sequence names identifying outliers.
    """

    assert aligner in {'cmalign', 'muscle', 'vsearch'}, 'invalid aligner: ' + aligner
    assert strategy in {'radius', 'cluster'}, 'invalid strategy: ' + strategy
    assert cluster_type in {'single', 'RobustSingleLinkage'}, \
        'invalid cluster_type ' + cluster_type

    prefix = '{}_'.format(tax_id)

    if distmat is None:
        log.debug('running {} on {}'.format(aligner, tax_id))
        if aligner == 'cmalign':
            taxa, distmat = distmat_cmalign(
                sequence_file, prefix, cpu=threads or wrap.CMALIGN_THREADS)
        elif aligner == 'muscle':
            taxa, distmat = distmat_muscle(sequence_file, prefix, maxiters)
        elif aligner == 'vsearch':
            taxa, distmat = distmat_pairwise(
                sequence_file, prefix, aligner, executable, iddef,
                threads=threads or wrap.VSEARCH_THREADS)
    else:
        assert taxa is not None

    if cutoff is not None:
        cutoff = cutoff
    elif percentile is not None:
        cutoff = outliers.scaled_radius(
            distmat, percentile, min_radius, max_radius)
    else:
        raise ValueError('must provide either cutoff or percentile')

    log.info('cutoff={} ({})'.format(
        cutoff, 'calculated' if percentile else 'pre-defined'))

    log.info('strategy: {}'.format(strategy))
    if strategy == 'radius':
        medoid, dists, is_out = outliers.outliers(distmat, cutoff)
        clusters = numpy.repeat(medoid, len(taxa))
    elif strategy == 'cluster':
        medoid, dists, is_out, clusters = outliers.outliers_by_cluster(
            distmat, t=cutoff, D=1.5,
            min_size=2, cluster_type=cluster_type)

    assert len(is_out) == len(taxa)

    result = pd.DataFrame({
        'seqname': taxa,
        'centroid': numpy.repeat(taxa[medoid], len(taxa)),
        'dist': dists,
        'is_out': is_out,
        'cluster': clusters})

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
                  fa_idx,
                  seqs,
                  strategy,
                  cluster_type,
                  distance_cutoff,
                  percentile,
                  min_radius,
                  max_radius,
                  aligner,
                  executable,
                  threads):
    """
    Worker task for running filtering tasks.

    Arguments:
    :tax_id: string identifying tax_id of these sequences
    :sequence_file: Complete sequence file
    :seqs: sequence names representing this tax_id
    :strategy: 'radius' or 'cluster'
    :distance_cutoff: Distance cutoff for medoid filtering
    :percentile: Used to calculate distance cutoff (``distance_cutoff`` cannot
     also be provided)
    :{min, max}_radius: Bounds of calculated distance cutoff
    :cluster_type: clustering algorithm (method of scipy.cluster.hierarchy)
    :aligner: name of alignment program
    :executable: name or path of executable for alignmment program

    :returns: output of ``filter_sequences()``
    """

    prefix = '{}_'.format(tax_id)

    with util.ntf(prefix=prefix, suffix='.fasta') as tf:
        # Extract sequences
        wrap.esl_sfetch(sequence_file, seqs, tf, fa_idx)
        tf.flush()

        filtered = filter_sequences(
            tax_id,
            sequence_file=tf.name,
            strategy=strategy,
            cluster_type=cluster_type,
            cutoff=distance_cutoff,
            percentile=percentile,
            min_radius=min_radius,
            max_radius=max_radius,
            aligner=aligner,
            executable=executable,
            threads=threads)

        return filtered


def action(a):
    # itemize sequences provided in the input file
    fa_idx = wrap.read_seq_file(a.sequence_file)
    seqnames = list(fa_idx.keys())

    # Load taxonomy
    with a.taxonomy as fp:
        taxonomy = TaxNode.from_taxtable(fp)
        log.info('Loaded taxonomy')

    # Load sequences into taxonomy
    with open(a.seqinfo_file) as fp:
        taxonomy.populate_from_seqinfo(fp, seqnames)
        n_added = sum(1 for i in taxonomy.subtree_sequence_ids())
        if n_added == 0:
            log.error('No sequences were added. Are all '
                      'tax_ids present in the taxonomy?')
            sys.exit(1)
        log.info('Added %d sequences', n_added)

    executable = a.executable or {'cmalign': 'cmalign',
                                  'muscle': 'muscle',
                                  'vsearch': wrap.VSEARCH}[a.aligner]

    if (a.previous_details and
            os.path.isfile(a.previous_details) and
            os.stat(a.previous_details).st_size):
        dtype = {'seqname': str, 'tax_id': str, a.filter_rank: str}
        # columns in output of `filter_worker`
        filter_worker_cols = [
            'centroid', 'cluster', 'dist', 'is_out', 'seqname', 'x', 'y']
        previous_details = pd.read_csv(
            a.previous_details, dtype=dtype).groupby(a.filter_rank)
    else:
        previous_details = None

    outcomes = []  # accumulate DatFrame objects

    # Sequences which are classified above the desired rank should just be kept
    names_above_rank = set(sequences_above_rank(taxonomy, a.filter_rank))
    log.info('Keeping %d sequences classified above %s',
             len(names_above_rank), a.filter_rank)
    above_rank = mock_filter(names_above_rank, keep=True)
    above_rank[a.filter_rank] = pd.Series(
        numpy.nan, index=above_rank.index)
    outcomes.append(above_rank)

    # For each filter-rank, filter
    nodes = [i for i in taxonomy if i.rank == a.filter_rank]

    # make sure all seqs are accounted for
    names_at_rank = set(s for n in nodes for s in n.subtree_sequence_ids())

    for s in seqnames:
        if s not in names_above_rank and s not in names_at_rank:
            raise ValueError(s + ' missing tax_id at filter rank')

    # --no-filter taxids
    if not a.no_filter:
        no_filter = set()
    elif os.path.isfile(a.no_filter):
        no_filter = set(i.strip() for i in open(a.no_filter) if i)
    else:
        no_filter = set(a.no_filter.split(','))

    # Filter each tax_id, running ``--jobs`` tasks in parallel
    with futures.ThreadPoolExecutor(a.jobs) as executor:
        # dispatch a pool of tasks
        futs = {}
        for i, node in enumerate(nodes):
            seqs = frozenset(node.subtree_sequence_ids())

            if previous_details and node.tax_id in previous_details.groups:
                prev_seqs = previous_details.get_group(node.tax_id)
            else:
                prev_seqs = None

            # in each case, `f` is a Future returning a DataFrame (see
            # filter_sequences)
            if not seqs:
                log.debug("No sequences for %s (%s)", node.tax_id, node.name)
                continue
            elif len(seqs) < a.min_seqs_for_filtering:
                log.debug('{} sequence(s) for {} ({}) [action: {}]'.format(
                    len(seqs), node.tax_id, node.name, a.rare_taxon_action))
                f = executor.submit(mock_filter, seqs=list(seqs),
                                    keep=a.rare_taxon_action == KEEP)
            elif node.tax_id in no_filter:
                log.debug('{} --no-filter sequence(s) for {} ({})'.format(
                    len(seqs), node.tax_id, node.name))
                f = executor.submit(mock_filter, seqs=list(seqs), keep=True)
            elif prev_seqs is not None and set(prev_seqs['seqname']) == seqs:
                # use previous results
                log.info(
                    'using previous results for tax_id {}'.format(node))
                f = executor.submit(
                    lambda f: f, f=prev_seqs[filter_worker_cols])
            else:
                f = executor.submit(
                    filter_worker,
                    tax_id=node.tax_id,
                    sequence_file=a.sequence_file,
                    fa_idx=fa_idx,
                    seqs=seqs,
                    strategy=a.strategy,
                    distance_cutoff=a.distance_cutoff,
                    percentile=a.distance_percentile,
                    min_radius=a.min_distance,
                    max_radius=a.max_distance,
                    cluster_type=a.cluster_type,
                    aligner=a.aligner,
                    executable=executable,
                    threads=a.threads_per_job)

            futs[f] = {'n_seqs': len(seqs), 'node': node}

        # log results for each tax_id as tasks complete
        complete = 0
        while futs:
            done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
            complete += len(done)
            sys.stderr.write('{0:8d}/{1:8d} taxa completed\r'.format(
                complete, complete + len(pending)))
            for f in done:
                exception = f.exception()
                if exception:
                    log.exception(
                        "Error in child process: %s", exception)
                    executor.shutdown(wait=False)
                    traceback.print_tb(f._traceback)
                    raise exception

                info = futs.pop(f)
                filtered = f.result().copy()  # here's the DataFrame again...

                # add a column for tax_d at filter_rank
                filtered.loc[:, a.filter_rank] = info['node'].tax_id
                outcomes.append(filtered)

                kept = frozenset(filtered.seqname[~filtered.is_out])
                if len(kept) == 0:
                    log.info('Pruned all %d sequences for %s (%s)',
                             info['n_seqs'], info['node'].tax_id,
                             info['node'].name)
                elif len(kept) != info['n_seqs']:
                    log.info('Pruned %d/%d sequences for %s (%s)',
                             info['n_seqs'] - len(kept), info['n_seqs'],
                             info['node'].tax_id, info['node'].name)

    all_outcomes = pd.concat(outcomes, ignore_index=True)
    all_outcomes.set_index('seqname', inplace=True)

    # all input sequences should be in the output
    assert {
        s for node in taxonomy for s in node.sequence_ids} == set(
        all_outcomes.index)

    kept_ids = set(all_outcomes.index[~all_outcomes.is_out])

    with a.output_seqs as fp:
        # Extract all of the sequences that passed.
        log.info('Extracting %d sequences', len(kept_ids))
        wrap.esl_sfetch(a.sequence_file, kept_ids, fp, fa_idx)

    # Filter seqinfo for sequences that passed.
    seqinfo = pd.read_csv(
        a.seqinfo_file, dtype={'seqname': str, 'tax_id': str})
    seqinfo = seqinfo.loc[seqinfo['seqname'].isin(seqnames)]
    seqinfo.set_index('seqname', inplace=True)

    merged = seqinfo.join(all_outcomes, lsuffix='.left')

    # csv output
    if a.filtered_seqinfo:
        merged[~merged.is_out].to_csv(
            a.filtered_seqinfo,
            columns=seqinfo.columns)

    if a.detailed_seqinfo:
        with open(a.detailed_seqinfo, 'w') as detailed_seqinfo:
            merged.to_csv(detailed_seqinfo)
