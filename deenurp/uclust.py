"""Wrapper and parsers for vsearch as a replacement for UCLUST.

Note that memory use can be high: searching against 2.2 million 16S
sequences from RDP takes ~30GB of memory for one strand search

TODO: rename this module and functions herein - keeping ``uclust`` for
now to avoid refactoring.

"""

import collections
import contextlib
import csv
import itertools
import logging
import operator
import subprocess
import tempfile

import numpy as np
import pandas as pd
from Bio import SeqIO

from .util import require_executable

log = logging.getLogger(__name__)

DEFAULT_PCT_ID = 0.99

# For parsing .uc format
UCLUST_HEADERS = ['type', 'cluster_number', 'size', 'pct_id',
                  'strand', 'query_start', 'seed_start', 'alignment',
                  'query_label', 'target_label']
UCLUST_TYPES = {'cluster_number': int, 'pct_id': float, 'query_start': int,
                'seed_start': int, 'size': int}

UClustRecord = collections.namedtuple('UClustRecord', UCLUST_HEADERS)


@contextlib.contextmanager
def _handle(s, *args, **kwargs):
    """
    Generate a file-like object from s.

    If s is a string, opens s and yields the open file. Otherwise, has no
    effect.
    """
    if isinstance(s, str):
        with open(s, *args, **kwargs) as fp:
            yield fp
    else:
        raise ValueError('try passing in a string instead')
        yield s


@contextlib.contextmanager
def _maybe_tempfile_name(obj=None, **kwargs):
    if obj:
        yield obj
    else:
        with tempfile.NamedTemporaryFile(**kwargs) as tf:
            yield tf.name


def _check_call(cmd, **kwargs):
    """
    Log and run command. Additional arguments are passed to
    ``subprocess.check_call``
    """
    cmd = list(map(str, cmd))
    logging.debug(' '.join(cmd))
    subprocess.check_call(cmd, **kwargs)


# Parsing
def _parse_uclust_row(row):
    if not len(row) == len(UCLUST_HEADERS):
        raise ValueError("Unexpected row length: {0} ({1})".format(len(row), row))
    for i, (header, val) in enumerate(zip(UCLUST_HEADERS, row)):
        # Replace NA char (*) with None, type convert
        if val in ('*', ''):
            row[i] = None
        elif header in UCLUST_TYPES:
            row[i] = UCLUST_TYPES[header](val)

    return UClustRecord(*row)


def parse_uclust_out(ucout_fp):
    """
    Parse the results of running UCLUST, returning UClustRecords.

    ucout_fp can be file name or text mode file handle.
    """

    with _handle(ucout_fp) as fp:
        # Skip comments
        rows = (i.rstrip() for i in fp if not i.startswith('#'))
        reader = csv.reader(rows, delimiter='\t')
        for row in reader:
            yield _parse_uclust_row(row)


def parse_uclust_as_df(ucout_fp):
    dtype = {'type': str, 'query_label': str,
             'target_label': str, 'alignment': str}
    df = pd.read_csv(
        ucout_fp, sep='\t', na_values='*', names=UCLUST_HEADERS, dtype=dtype)

    # define target_label as query_label for seed sequences
    df['target_label'] = np.where(
        df['type'] == 'S', df['query_label'], df['target_label'])

    return df


# Library search
def hits_by_sequence(uclust_records):
    """
    Collect hits by sequence

    Generates (sequence_id, list_of_hits) tuples from an iterable of uclust
    records.

    uclust_records must be sorted such that hits for a single sequences occur
    sequentially (uclust output fulfills this without modification).
    """
    # Only interested in hit (H) and NoHit (N) records
    uclust_records = (i for i in uclust_records if i.type in ('H', 'N'))
    grouped = itertools.groupby(uclust_records, operator.attrgetter('query_label'))
    for g, v in grouped:
        # Only return hits
        yield g, [i for i in v if i.type == 'H']


def sequences_by_cluster(uclust_records):
    """
    Collect sequences by cluster

    Generates (seed_sequence_id, records_in_cluster) tuples from an iterable
    of uclust records.
    """
    # Seeds and hits
    k = operator.attrgetter('cluster_number')
    r = (i for i in uclust_records if i.type in ('S', 'H'))
    s = sorted(r, key=k)
    grouped = itertools.groupby(s, k)
    for g, v in grouped:
        l = list(v)
        yield l[0].query_label, l


def cluster_map(uclust_records):
    """
    Map sequence names to clusters.

    Generates (cluster_number, sequence_name, seed_sequence_name) tuples
    from an iterable of uclust records.
    """

    for row in uclust_records:
        if row.type == 'S':
            yield (row.cluster_number, row.query_label, row.query_label)
        elif row.type == 'H':
            yield (row.cluster_number, row.query_label, row.target_label)


def search(database, query, output, pct_id=DEFAULT_PCT_ID,
           maxaccepts=None, maxrejects=None, quiet=False, search_pct_id=None):
    """
    Run UCLUST against a sequence database in FASTA format.

    Parameters:
     database:        Path to FASTA file to search against
     query:           Path to query file
     output:          Path for uclust output
     pct_id:          Minimum identity for match (provided to ``uclust --id``)
     search_pct_id:   If given, the database is searched at search_pct_id, then
                      the results filtered to only include sequences that match
                      at greater than pct_id.

                      Note: If search_pct_id is specified, cluster sizes will
                      be inaccurate.

    Others: see ``vsearch --help``
    """
    require_executable('vsearch')
    with _maybe_tempfile_name(
            output if not search_pct_id else None, prefix='vsearch-') as o:
        cmd = ['vsearch',
               '--usearch_global', query,
               '--db', database,
               '--uc', o,
               '--uc_allhits',  # show all, not just top hit with uc output
               '--id', str(search_pct_id or pct_id)]  # Prefer search_pct_id
        if maxaccepts:
            cmd.extend(('--maxaccepts', str(maxaccepts)))
        if maxrejects:
            cmd.extend(('--maxrejects', str(maxrejects)))
        if quiet:
            cmd.append('--quiet')

        _check_call(cmd)

        if search_pct_id:
            # Filter results, write to output
            with open(output, 'w') as uc:
                w = csv.writer(uc, lineterminator='\n', delimiter='\t')
                records = parse_uclust_out(o)
                # Filter records which don't meet the pct_id criteria
                # IDs in the output file are reported as percentages.
                id_cutoff = pct_id * 100.0
                records = (i for i in records if i.type != 'H' or i.pct_id >= id_cutoff)
                w.writerows(records)


def cluster(sequence_file, output, pct_id=DEFAULT_PCT_ID, quiet=False,
            pre_sorted=False, threads=None):
    """Cluster de novo. If ``pre_sorted`` is True, assume that sequences
    are pre-sorted by length (and cluster using --cluster_smallmem
    rather than --cluster_fast). See ``vsearch --help`` for details.

    """
    require_executable('vsearch')
    cmd = ['vsearch',
           '--cluster_smallmem' if pre_sorted else '--cluster_fast', sequence_file,
           '--uc', output,
           '--id', str(pct_id)]
    if quiet:
        cmd.append('--quiet')
    if not pre_sorted:
        cmd.append('--usersort')
    if threads is not None:
        cmd.extend(['--threads', str(threads)])
    _check_call(cmd)


def cluster_seeds(sequence_file, uclust_out):
    """
    Extract seeds from the result of clustering, e.g.:

    >>> cluster('seqs.fasta', 'seqs.99.uc')
    >>> seeds = cluster_seeds('seqs.fasta', 'seqs.99.uc')
    """
    uc = parse_uclust_out(uclust_out)
    seeds = frozenset(i.query_label.split(None, 1)[0] for i in uc if i.type == 'S')
    seen_seeds = set()

    sequences = SeqIO.parse(sequence_file, 'fasta')

    for s in sequences:
        if s.id in seeds:
            seen_seeds.add(s.id)
            yield s

    if seeds - seen_seeds:
        raise ValueError(
            "Some expected seeds were not found in the FASTA file: {0}".format(
                ','.join(seeds - seen_seeds)))

# Functions to convert uclust output into format usable by `guppy redup -m`


class DeduplicatedSequence(object):
    __slots__ = ['id', 'count']

    def __init__(self, id, count):
        self.id = id
        self.count = count


class ConstantDict(object):
    def __init__(self, return_value=None):
        self.return_value = return_value

    def __getitem__(self, item):
        return self.return_value


def guppy_redup_from_uclust(uclust_records, sample_map=None):
    """
    Generate a guppy-redup compatible mapping file from UCLUST output.

    If sample map is specified, at one sequence is kept in the output map from
    every specimen present in each cluster, along with a count of sequences
    from the specimen within the cluster.

    Otherwise, each cluster is reduced to the seed id and a count of sequences
    within.
    """
    sample_map = sample_map or ConstantDict()
    # Map from cluster number to seed id
    seeds = {}
    # Map from cluster number to map of sample->seed,count
    clusters = collections.defaultdict(dict)
    # Just seeds and hits
    uclust_records = (i for i in uclust_records if i.type in ('S', 'H'))

    for record in uclust_records:
        number = record.cluster_number
        q = record.query_label.split()[0]
        sample = sample_map[q]
        if record.type == 'S':
            seeds[number] = q
        if sample not in clusters[number]:
            clusters[number][sample] = DeduplicatedSequence(q, 1)
        else:
            clusters[number][sample].count += 1

    rows = [(seeds[num], dedup_seq.id, dedup_seq.count)
            for num, samples in list(clusters.items())
            for dedup_seq in list(samples.values())]
    return rows
