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

    # TODO: there has to be a more efficient way to read a numpy matrix
    elements = fobj.read().split()
    N = int(elements.pop(0))
    distmat = numpy.repeat(-1*numpy.inf, N**2)
    distmat.shape = (N, N)

    taxa = []
    for row, i in enumerate(range(0, len(elements), N+1)):
        taxa.append(elements[i])
        distmat[row, :] = [float(x) for x in elements[i+1:i+N+1]]

    return taxa, distmat

def fasttree_dists(fasta, suppress_stderr = True):
    """
    Calculate pairwise distances among DNA multiple alignment in
    `fasta` using FastTree and return (taxon_names, matrix).
    """

    # TODO: need a more informative error if FastTree is not installed.

    cmd = ['FastTree','-nt','-makematrix', fasta]

    with tempfile.TemporaryFile('rw') as stdout, open(os.devnull) as devnull:
        proc = subprocess.Popen(
            cmd, stdout = stdout,
            stderr = devnull if suppress_stderr else None)
        proc.communicate()
        stdout.flush()
        stdout.seek(0)
        taxa, distmat = read_dists(stdout)

    return taxa, distmat

def outliers(mat, cutoff):
    """
    Given pairwise distance matrix `mat`, identify elements with a
    distance to the centrid element of > cutoff. Returns a boolean
    vector corresponding to the margin of mat.
    """

    # index of most central element.
    medoid = numpy.argmin(numpy.median(mat, 0))

    # distance from each element to most central element
    dists = mat[medoid, :]

    return dists > cutoff
