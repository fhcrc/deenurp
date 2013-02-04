"""
Identify mislabeled sequence records.
"""

import os
import logging
import subprocess
import tempfile

import numpy

log = logging.getLogger(__name__)

def read_dists(fobj):
    """
    Read interleaved phylip distance matrix from file-like object fobj
    into a numpy matrix. Return (taxon_names, matrix).
    """

    N = int(fobj.readline())
    distmat = numpy.repeat(-1*numpy.inf, N**2)
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

    cmd = ['FastTree','-nt','-makematrix', fasta]

    with tempfile.TemporaryFile('rw') as stdout, open(os.devnull) as devnull:
        proc = subprocess.Popen(cmd, stdout = stdout, stderr = devnull)
        proc.communicate()
        stdout.flush()
        stdout.seek(0)
        taxa, distmat = read_dists(stdout)

    return taxa, distmat

def outliers(distmat, cutoff):
    """
    Given pairwise distance matrix `mat`, identify elements with a
    distance to the centrid element of > cutoff. Returns a boolean
    vector corresponding to the margin of mat.
    """

    if distmat.shape[0] == 2:
        # it really doesn't make sense to perform this test at all
        # with a small number of comparisons, but we'll try to make
        # some sense out of a comparison of only two sequences by
        # discarding both if the distance between them is > cutoff.
        both_fail = distmat[0, 1] > cutoff
        return numpy.repeat(both_fail, 2)
    else:
        # index of most central element.
        medoid = numpy.argmin(numpy.median(distmat, 0))

        # distance from each element to most central element
        dists = distmat[medoid, :]

        return dists > cutoff
