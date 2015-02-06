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

log = logging.getLogger(__name__)


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
        distmat[row, :] = map(float, spl)

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
    """Return the index of the medoid of square numpy matrix
    ``X``. ``ii`` is an optional boolean vector specifying which rows
    and columns of ``X`` to consider.

    """

    n, m = X.shape
    assert n == m

    if ii is None:
        medoid = np.argmin(np.median(X, 0))
        idx = np.arange(n)[medoid]
    else:
        medoid = np.argmin(np.median(X[np.ix_(ii, ii)], 0))
        idx = np.arange(n)[ii][medoid]

    return idx


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


def outliers_by_cluster(distmat, t, max_dist, min_size=2, cluster_type='single'):

    clusters, title = scipy_cluster(distmat, cluster_type, t=t)
    medoids = find_cluster_medoids(distmat, clusters)
    keep = choose_clusters(medoids, min_size, max_dist)
    to_prune = ~ pd.Series(clusters).isin(keep)

    # medoid of the largest cluster
    medoid = medoids['medoid'][0]
    dists = None

    return medoid, dists, to_prune


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
        module, ' '.join('%s=%s' % item for item in args.items()))

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
    assert n == m, 'X must be a square matrix'
    assert isinstance(clusters, np.ndarray)

    clusters, counts = np.unique(clusters, return_counts=True)
    tallies = sorted(zip(counts, clusters), reverse=True)
    counts, clusters = zip(*tallies)  # reorders clusters and counts
    medoids = [(None if cluster == -1 else find_medoid(X, clusters == cluster))
               for _, cluster in tallies]
    dists = [None if medoid is None else X[medoids[0], medoid] for medoid in medoids]

    return pd.DataFrame.from_items([
        ('cluster', clusters), ('count', counts), ('medoid', medoids), ('dist', dists)
    ]).sort('count', ascending=False)


def choose_clusters(df, min_size, max_dist):
    """Implements logic for cluster-based outlier detection.

    * ``df`` - output of find_cluster_medoids()
    * ``min_size`` - discard clusters with size below this value.
    * ``max_dist`` - discard clusters whose medoid is greater than
      this distance from the medoid of the largest cluster.

    In addition to selecting clusters according to the parameters
    above, also discards any clusters identified by a value of -1.

    Returns an ndarray containing names clusters to keep (ie, values
    from 'cluster' column).

    """

    # the parens are necessary to enforce precedence
    keep = (df['cluster'] != -1) & (df['count'] >= min_size) & (df['dist'] <= max_dist)
    return df['cluster'][keep]


def scaled_radius(X, percentile, min_radius=0.0):
    """Calculate the distribution of distances between the medoid of
    distance matrix ``X`` and each element, and return the value at
    ``percentile``. ``min_radius`` defines a minimum return value.

    """

    return max([np.percentile(X[find_medoid(X), :], percentile), min_radius])
