"""
Identify mislabeled sequence records.
"""

import os
import logging
import subprocess
import tempfile

import numpy as np

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


def outliers(distmat, cutoff, min_size=3):
    """Given pairwise distance matrix `distmat`, identify elements with a
    distance to the centrid element of > cutoff. Does not attempt to
    prune if margin of distmat is < min_size.

    Returns (medoid, dists, to_prune):

    * medoid - the index of the centermost element
    * dists - a vector of distances to the medoid
    * to_prune - a boolean vector corresponding to the margin of `distmat`
      where True identifies outliers.

    """

    if distmat.shape[0] < min_size:
        medoid = np.nan
        dists = np.repeat(np.nan, distmat.shape[0])
        to_prune = np.repeat(False, distmat.shape[0])
    else:
        # use a masked array in case there are any nan
        ma = np.ma.masked_array(distmat, np.isnan(distmat))

        # index of most central element.
        medoid = np.argmin(np.median(ma, 0))

        # distance from each element to most central element
        dists = ma[medoid, :]
        to_prune = dists > cutoff

    return medoid, dists, to_prune
