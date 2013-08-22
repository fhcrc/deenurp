"""
Wrapper and parsers for UCLUST.

UCLUST is a precursor to USEARCH, released freely with PyNAST and QIIME. The
binary we have is 64-bit, so large database searches are possible.

Note that memory use can be high: searching against 2.2 million 16S
sequences from RDP takes ~30GB of memory for one strand search
"""

import os.path
import collections
import contextlib
import csv
import itertools
import logging
import operator
import subprocess
import tempfile

from Bio import SeqIO

from .util import require_executable

log = logging.getLogger(__name__)

DEFAULT_PCT_ID = 0.99

# For parsing .uc format
UCLUST_HEADERS = ['type', 'cluster_number', 'size', 'pct_id', 'strand',
    'query_start', 'seed_start', 'alignment', 'query_label', 'target_label']
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
    if isinstance(s, basestring):
        with open(s, *args, **kwargs) as fp:
            yield fp
    else:
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
    cmd = map(str, cmd)
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

    ucout_fp can be file name or file handle.
    """
    with _handle(ucout_fp) as fp:
        # Skip comments
        rows = (i.rstrip() for i in fp if not i.startswith('#'))
        reader = csv.reader(rows, delimiter='\t')
        for row in reader:
            yield _parse_uclust_row(row)


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
    grouped = itertools.groupby(uclust_records,
            operator.attrgetter('query_label'))
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
        maxaccepts=None, maxrejects=None, quiet=False, trunclabels=False,
        wordcountreject=True, search_pct_id=None):
    """
    Run UCLUST against a sequence database in FASTA format.

    Parameters:
     database:        Path to FASTA file to search against
     query:           Path to query file
     output:          Path for uclust output
     pct_id:          Minimum percent ID for match
     wordcountreject: Pre-filter based on word count in common between query
                      and target seqs. When true, decreases sensitivity, but
                      decreases speed by 50-70%.
     search_pct_id:   If given, the database is searched at search_pct_id, then
                      the results filtered to only include sequences that match
                      at greater than pct_id. *May* produce the same results as
                      disabling wordcountreject, will be faster.

                      Note: If search_pct_id is specified, cluster sizes will
                      be inaccurate.

    Others: see ``uclust --help``
    """
    require_executable('uclust')
    with _maybe_tempfile_name(output if not search_pct_id else None, prefix='uclust-') as o:
        cmd = ['uclust',
                '--input', query,
                '--lib', database,
                '--uc', o,
                '--libonly',       # Don't generate new clusters when no DB hits
                '--allhits',       # output all hits
                '--id', str(search_pct_id or pct_id)]  # Prefer search_pct_id
        if not wordcountreject:
            cmd.append('--nowordcountreject')
        if maxaccepts:
            cmd.extend(('--maxaccepts', str(maxaccepts)))
        if maxrejects:
            cmd.extend(('--maxrejects', str(maxrejects)))
        if quiet:
            cmd.append('--quiet')
        if trunclabels:
            cmd.append('--trunclabels')

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
        usersort=False, trunclabels=False, wordcountreject=True):
    """
    Cluster de novo
    """
    require_executable('uclust')
    cmd = ['uclust',
            '--input', sequence_file,
            '--uc', output,
            '--id', str(pct_id)]
    if not wordcountreject:
        cmd.append('--nowordcountreject')
    if quiet:
        cmd.append('--quiet')
    if usersort:
        cmd.append('--usersort')
    if trunclabels:
        cmd.append('--trunclabels')
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
            "Some expected seeds were not found in the FASTA file: {0}".format(','.join(seeds - seen_seeds)))

def sort(sequence_file, output, quiet=False):
    """
    sort by descending length
    """
    require_executable('uclust')
    cmd = ['uclust', '--sort', sequence_file, '--output', output]
    if quiet:
        cmd.append('--quiet')
    _check_call(cmd)


def sort_and_cluster(sequence_file, output, **kwargs):
    """
    Sort sequence_file by descending length, then cluster.

    Additional arguments are passed to cluster
    """
    with tempfile.NamedTemporaryFile(prefix=os.path.basename(sequence_file)) as tf:
        sort(sequence_file, tf.name, kwargs.get('quiet'))
        cluster(tf.name, output, **kwargs)

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
        if not sample in clusters[number]:
            clusters[number][sample] = DeduplicatedSequence(q, 1)
        else:
            clusters[number][sample].count += 1

    rows = [(seeds[number], dedup_seq.id, dedup_seq.count)
            for number, samples in clusters.items()
            for dedup_seq in samples.values()]
    return rows
