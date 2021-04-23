"""
Identify mislabeled sequence records.
"""

import os
import logging
import subprocess
import tempfile

import numpy as np
import pandas as pd

import scipy
import scipy.cluster

import hdbscan

log = logging


def read_dists(fobj):
    """
    Read interleaved phylip distance matrix from file-like object fobj
    into a numpy matrix. Return (taxon_names, matrix).
    """

    N = int(fobj.readline())
    distmat = np.repeat(-1 * np.inf, N ** 2)
    distmat.shape = (N, N)

    taxa = []
    for row, line in enumerate(fobj):
        spl = line.split()
        assert len(spl) == N + 1
        taxa.append(spl.pop(0))
        distmat[row, :] = list(map(float, spl))

    return taxa, distmat


def fasttree_dists(fasta):
    """
    Calculate pairwise distances among DNA multiple alignment in
    `fasta` using FastTree and return (taxon_names, matrix).
    """

    # TODO: need a more informative error if FastTree is not installed.

    cmd = ['FastTree', '-nt', '-makematrix', fasta]

    with tempfile.TemporaryFile('rw') as stdout, open(os.devnull) as devnull:
        proc = subprocess.Popen(cmd, stdout=stdout, stderr=devnull)
        proc.communicate()
        stdout.flush()
        stdout.seek(0)
        taxa, distmat = read_dists(stdout)

    return taxa, distmat


def find_medoid(X, ii=None):
    """Return the index of the medoid of square numpy matrix ``X`` of
    shape (n, n). ``ii`` is an optional boolean vector of length n
    defining a submatrix of ``X``.

    """

    n, m = X.shape
    assert n == m, 'X must be a square matrix'

    if ii is None:
        medoid = np.argmin(np.median(X, 0))
    else:
        assert ii.shape[0] == n, 'ii must be the length of the margin of X'
        medoid_ii = np.argmin(np.median(X[np.ix_(ii, ii)], 0))
        medoid = np.arange(n)[ii][medoid_ii]

    return medoid


def all_ok(distmat):
    medoid = np.nan
    dists = np.repeat(np.nan, distmat.shape[0])
    to_prune = np.repeat(False, distmat.shape[0])
    return medoid, dists, to_prune


def outliers(distmat, radius):
    """Given pairwise distance matrix `distmat`, identify elements where
    distance to the medoid > radius.

    Returns (medoid, dists, to_prune):

    * medoid - the index of the centermost element
    * dists - a vector of distances to the medoid
    * to_prune - a boolean vector corresponding to the margin of `distmat`
      where True identifies outliers.

    """

    # use a masked array in case there are any nan
    ma = np.ma.masked_array(distmat, np.isnan(distmat))

    # index of most central element.
    medoid = find_medoid(ma)

    # distance from each element to most central element
    dists = ma[medoid, :]
    to_prune = dists > radius

    return medoid, dists, to_prune


def outliers_by_cluster(distmat, t, D, min_size=1, cluster_type='single', **kwargs):
    """Detect outliers by 1) performing hierarchical clustering based on
    distances in ``distmat`` (a square numpy matrix) with distance
    threshold ``t`` and discarding clusters with fewer than
    ``min_size`` members; then 2) discarding clusters whose medoids
    have a distance greater than ``t * D`` from the medoid of the
    largest cluster. ``cluster_type`` is a string naming a submodule
    of either ``scipy.cluster.hierarchy`` or ``hdbscan`` (ie,
    'HDBSCAN' or 'RobustSingleLinkage') and identifies the clustering
    algorithm. Additional parameters are passed to the clustering
    algorithm when using one of the hdbscan methods only.

    Return values:

    * medoid - the index of the centermost element
    * dists - a vector of distances to the medoid of the largest cluster
    * to_prune - a boolean vector corresponding to the margin of `distmat`
      where True identifies outliers.
    * clusters - an index into `distmat` representing each cluster

    """

    if cluster_type == 'HDBSCAN':
        clusters, title = hdbscan_cluster(distmat, cluster_type, **kwargs)
    elif cluster_type == 'RobustSingleLinkage':
        clusters, title = hdbscan_cluster(distmat, cluster_type, cut=t, **kwargs)
    else:
        clusters, title = scipy_cluster(distmat, cluster_type, t=t)

    log.info(title)

    if all(clusters == -1):
        # if no clusters are found, we throw up our hands and accept
        # all of the sequences
        log.warning('no clusters were found')

        medoids = pd.DataFrame.from_items([
            ('cluster', [-1]),
            ('count', [len(clusters)]),
            ('medoid', [find_medoid(distmat)]),
            ('dist', [None])
        ])
        to_prune = pd.Series([False for x in clusters])
    else:
        medoids = find_cluster_medoids(distmat, clusters)

        # `keep` is a collection of cluster numbers to retain (these are
        # cluster labels, not indices)
        keep = choose_clusters(medoids, min_size, t * D)
        to_prune = ~ pd.Series(clusters).isin(keep)

    # medoid of the largest cluster (an index into distmat)
    medoid = int(medoids['medoid'][0])

    # distances to this medoid (note that these distances are not used
    # directly as a criterion for filtering - they just provide a
    # sense of the relative position of each sequence)
    dists = distmat[medoid, :]

    return medoid, dists, to_prune, clusters


def scipy_cluster(X, module, t, **kwargs):
    """Given distance matrix ``X``, perform clustering using the specified
    module of ``scipy.cluster.hierarchy``. The threshold value `t` and
    any additional keyword arguments are passed to
    scipy.cluster.hierarchy.fcluster().

    see https://docs.scipy.org/doc/scipy-0.15.1/reference/cluster.hierarchy.html

    """

    defaults = {'criterion': 'distance'}
    args = dict(defaults, **kwargs)

    # requires 'import scipy.cluster'
    fun = getattr(scipy.cluster.hierarchy, module)
    y = scipy.spatial.distance.squareform(X)
    Z = fun(y)
    clusters = scipy.cluster.hierarchy.fcluster(Z, t, **args)
    title = 'scipy.cluster.hierarchy.{} {}'.format(
        module, ' '.join('%s=%s' % item for item in list(args.items())))

    return clusters, title


def hdbscan_cluster(X, module, **kwargs):

    try:
        fun = getattr(hdbscan, module)
    except AttributeError:
        raise ValueError("'module' must be 'HDBSCAN' or 'RobustSingleLinkage'")

    clusterer = fun(metric='precomputed', **kwargs)

    title = ' '.join(str(clusterer).split())
    clusters = clusterer.fit_predict(X)

    return clusters, title


def find_cluster_medoids(X, clusters):
    """Inputs are ``X``, a square distance matrix, and ``clusters``, a
    1-dimensional array assigning each element in ``X`` to a
    cluster. Returns a pandas DataFrame with rows corresponding to
    clusters (sorted by size, descending) with columns 'cluster',
    'count', 'medoid' (an int identifying the medoid of each cluster),
    and 'dist' (distance from the medoid of this cluster to the medoid
    of the largest cluster).

    """

    assert isinstance(X, np.ndarray)
    n, m = X.shape
    assert n == m, '`X` must be a square matrix'
    assert isinstance(clusters, np.ndarray), '`clusters` must be a numpy ndarray'
    assert clusters.shape[0] == n, '`clusters` must be the size of the margin of `X`'

    uclusters, counts = np.unique(clusters, return_counts=True)

    # reorders uclusters and counts by decreasing size, placing tally
    # of unclustered elements last
    tallies = sorted(
        zip([0 if c == -1 else 1 for c in uclusters], counts, uclusters),
        reverse=True)

    __, counts, uclusters = list(zip(*tallies))

    medoids = [(None if cluster == -1 else find_medoid(X, clusters == cluster))
               for _, _, cluster in tallies]

    # measure distances from the medoid of the first (largest) cluster
    dists = [None if medoid is None else X[medoids[0], medoid] for medoid in medoids]

    return pd.DataFrame.from_items([
        ('cluster', uclusters),
        ('count', counts),
        ('medoid', medoids),
        ('dist', dists)
    ])


def choose_clusters(df, min_size, max_dist):
    """Implements logic for cluster-based outlier detection.

    Arguments:
    * ``df`` - output of find_cluster_medoids()
    * ``min_size`` - discard clusters with size below this value.
    * ``max_dist`` - discard clusters whose medoid is greater than
      this distance from the medoid of the largest cluster.

    In addition to selecting clusters according to the parameters
    above, also discards any clusters identified by a value of -1
    (identifies unassigned items in the output of some clustering
    methods).

    Returns an ndarray identifying clusters to keep (ie, values from
    'cluster' column).

    """

    log.info('\n' + str(df))

    # the parens are necessary to enforce precedence
    keep = (df['cluster'] != -1) & (df['count'] >= min_size) & (df['dist'] <= max_dist)
    return df['cluster'][keep]


def scaled_radius(X, percentile, min_radius=0.0, max_radius=None):
    """Calculate the distribution of distances between the medoid of
    distance matrix ``X`` and each element, and return the value at
    ``percentile``. ``min_radius`` and ``max_radius`` define a minimum
    and maximum return values, respectively.

    """

    radius = np.percentile(X[find_medoid(X), :], percentile)
    log.info('calculated cutoff: {}'.format(radius))
    if radius < min_radius:
        radius = min_radius

    if max_radius is not None and radius > max_radius:
        radius = max_radius

    log.info('cutoff after bounds check: {}'.format(radius))
    return radius


def mds(X, taxa, n_jobs=1):
    """Perform multidimensional scaling using ``sklearn.manifold`` given
    square distance matrix ``X``. Return a DataFrame with columns
    'seqname', 'x', 'y' in which 'seqname' contains the names provided
    in `taxa`.

    """

    n, m = X.shape
    assert n == m, 'X must be a square matrix'

    from sklearn import manifold
    mds = manifold.MDS(
        dissimilarity='precomputed',
        random_state=12345,
        n_jobs=n_jobs)

    if np.all(X == 0):
        df = pd.DataFrame.from_items([
            ('seqname', taxa),
            ('x', np.zeros(n)),
            ('y', np.zeros(n))
        ])
    else:
        mds_fit = mds.fit_transform(X)
        df = pd.DataFrame(mds_fit, columns=['x', 'y'])
        df['seqname'] = pd.Series(taxa, index=df.index)
        df = df[['seqname', 'x', 'y']]  # reorder columns

    assert X.shape[0] == df.shape[0]
    return df
