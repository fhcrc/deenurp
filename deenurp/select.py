"""
Select reference sequences for inclusion
"""
import collections
import csv
import functools
import itertools
import logging
import operator
import tempfile
import hashlib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import search, uclust
from concurrent import futures

from . import util, wrap
from .config import DEFAULT_THREADS
from .util import as_fasta, tempdir
from .wrap import (cmalign, as_refpkg, redupfile_of_seqs,
                   rppr_min_adcl, guppy_redup, pplacer, esl_sfetch)

CLUSTER_THRESHOLD = 0.998

"""
Minimum proportion of total mass in a cluster
to require before including references
"""
MIN_CLUSTER_PROP = -1.0


def log_error(fn):
    """
    Log errors but still raise them
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            logging.exception("ERROR in {}".format(fn.__name__))
            raise
    return wrapper


def seqrecord(name, residues, **annotations):
    """
    Takes a seqname and sequence plus annotations, constucts
    and return a Bio.SeqRecord.SeqRecord
    """
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr


def _cluster(sequences, threshold=CLUSTER_THRESHOLD):
    """
    Cluster ``sequences`` at ``threshold``, returning only the seeds.

    see CLUSTER_THRESHOLD
    """
    sequences = list(sequences)
    assert sequences
    with as_fasta(sequences) as fasta_name, \
            tempfile.NamedTemporaryFile(prefix='uc-') as ntf:

        uclust.cluster(fasta_name, ntf.name, pct_id=threshold, quiet=True)
        ntf.seek(0)
        r = list(uclust.cluster_seeds(fasta_name, ntf))

    logging.debug("Clustered %d to %d", len(sequences), len(r))
    return r


@log_error
def select_sequences_for_cluster(
        ref_seqs,
        query_seqs,
        cluster_name,
        cluster_weight,
        max_sample,
        max_weight,
        norm_sw,
        keep_leaves=5):
    """
    Given a set of reference sequences and query sequences, select
    keep_leaves appropriate references.
    """
    logging.info('Cluster %s: Max sample abundance: %.3f%% of %s, %d hits',
                 cluster_name, max_weight * 100, max_sample, len(query_seqs))

    # cluster to minimize redundancy among refs
    clustered = _cluster(ref_seqs, threshold=CLUSTER_THRESHOLD)

    # address the edge case in which only a single reference remains
    # after clustering at the above threshold:
    if len(clustered) < 3:
        clustered = _cluster(ref_seqs, threshold=1.0)

    ref_seqs = clustered

    mean_weight = sum(norm_sw.values()) / len(norm_sw)
    for ref in ref_seqs:
        ref.annotations.update({'cluster_name': cluster_name,
                                'max_weight': max_weight,
                                'mean_weight': mean_weight})

    if len(ref_seqs) <= keep_leaves:
        return ref_seqs
    # else: find some more reps

    # the operation below assumes unique identifiers for the set of
    # ref and query seqs, so ensure that this is the case
    for seq in query_seqs:
        seq.id = seq.id + hashlib.md5(seq.id).hexdigest()[:8]

    c = itertools.chain(ref_seqs, query_seqs)

    ref_ids = frozenset(i.id for i in ref_seqs)
    aligned = list(cmalign(c))
    with as_refpkg((i for i in aligned if i.id in ref_ids)) as rp, \
            as_fasta(aligned) as fasta, \
            tempdir(prefix='jplace') as placedir, \
            redupfile_of_seqs(query_seqs) as redup_path:

        jplace = pplacer(rp.path, fasta, out_dir=placedir(), threads=1)
        # Redup
        guppy_redup(jplace, redup_path, placedir('redup.jplace'))
        prune_leaves = set(
            rppr_min_adcl(placedir('redup.jplace'), keep_leaves))

    result = frozenset(i.id for i in ref_seqs) - prune_leaves
    assert len(result) == keep_leaves

    refs = [i for i in ref_seqs if i.id in result]
    return refs


@log_error
def select_sequences_for_whitelist_cluster(
        ref_seqs, cluster_name, keep_leaves=5):
    """
    Selects a subset of ``ref_seqs`` using ``rppr min_adcl_tree``
    """
    logging.info("Whitelisted cluster {}".format(cluster_name))

    # shrink ref_seqs by clustering first at 99.8% (CLUSTER_THRESHOLD)
    ref_seqs = _cluster(ref_seqs, threshold=CLUSTER_THRESHOLD)
    ref_ids = set(i.id for i in ref_seqs)
    for ref in ref_seqs:
        ref.annotations.update({'cluster_name': cluster_name,
                                'max_weight': None,
                                'mean_weight': None})

    if len(ref_seqs) <= keep_leaves:
        return ref_seqs

    aligned = list(cmalign(ref_seqs))
    with util.ntf(suffix='.tre') as tf:
        wrap.fasttree(aligned, tf, gtr=True)
        tf.close()
        prune = wrap.rppr_min_adcl_tree(tf.name, keep_leaves)

    result = ref_ids - frozenset(prune)
    assert len(result) == keep_leaves

    refs = [i for i in ref_seqs if i.id in result]
    return refs


def fetch_cluster_members(cluster_info_file, group_field):
    """
    """
    d = collections.defaultdict(list)
    with open(cluster_info_file) as fp:
        r = csv.DictReader(fp)
        for i in r:
            d[i[group_field]].append(i['seqname'])
    return dict(d)


def get_sample_weights(con, sequence_names):
    """
    Map from sample -> total weight for all samples associated with
    sequence_names
    """
    weights = collections.defaultdict(float)
    cursor = con.cursor()
    for sequence in sequence_names:
        sql = """SELECT samples.name AS sample_name, weight
        FROM sequences s
        INNER JOIN sequences_samples ss USING (sequence_id)
        INNER JOIN samples USING (sample_id)
        WHERE s.name = ?"""
        logging.debug(sql.replace('?', '{}').format(sequence))
        cursor.execute(sql, [sequence])
        res = cursor.fetchone()
        if not res:
            continue
        sample, weight = res
        weights[sample] += weight
    return dict(weights)


def sequences_hitting_cluster(con, cluster_name):
    """
    """
    sql = """SELECT DISTINCT sequences.name
        FROM sequences
        INNER JOIN best_hits USING (sequence_id)
        INNER JOIN ref_seqs USING(ref_id)
        WHERE cluster_name = ?
        ORDER BY sequences.name"""
    cursor = con.cursor()
    logging.debug(sql.replace('?', '{}').format(cluster_name))
    cursor.execute(sql, [cluster_name])
    return [name for name, in cursor]


def esl_sfetch_seqs(sequence_file, sequence_names, fa_idx):
    """
    """
    with tempfile.NamedTemporaryFile(prefix='esl', suffix='.fasta') as tf:
        esl_sfetch(sequence_file, sequence_names, tf, fa_idx)
        tf.seek(0)
        return list(SeqIO.parse(tf, 'fasta'))


def get_total_weight_per_sample(con):
    """
    """
    sql = """SELECT name, SUM(weight) AS weight
    FROM sequences_samples
      INNER JOIN samples USING (sample_id)
    GROUP BY name"""
    cursor = con.cursor()
    logging.debug(sql)
    cursor.execute(sql)
    return dict(cursor)


def choose_references(
        deenurp_db,
        ref_idx,
        fa_idx,
        refs_per_cluster=5,
        threads=DEFAULT_THREADS,
        min_cluster_prop=MIN_CLUSTER_PROP,
        include_clusters=None,
        exclude_clusters=None,
        include_sequences=None,
        exclude_sequences=None):
    """
    Choose reference sequences from a search, choosing refs_per_cluster
    reference sequences for each nonoverlapping cluster.

    min_cluster_prop - Minimum proportion of total mass in a cluster to
                       require before including references
    """

    if include_sequences:
        raise NotImplementedError('"include_sequences" is not implemented')

    params = search.load_params(deenurp_db)
    fasta_file = params['fasta_file']
    ref_fasta = params['ref_fasta']
    sample_total_weights = get_total_weight_per_sample(deenurp_db)
    cluster_members = fetch_cluster_members(
        params['ref_meta'], params['group_field'])

    # Iterate over clusters
    cursor = deenurp_db.cursor()

    # Select all clusters above cutoff
    sql = """
SELECT ref_seqs.cluster_name,
       sample_id,
       SUM(sequences_samples.weight) AS total_weight
FROM ref_seqs
    INNER JOIN best_hits USING (ref_id)
    INNER JOIN sequences USING (sequence_id)
    INNER JOIN sequences_samples USING (sequence_id)
GROUP BY ref_seqs.cluster_name, sample_id
HAVING SUM(sequences_samples.weight) > ?
ORDER BY ref_seqs.cluster_name ASC, SUM(sequences_samples.weight) DESC
"""
    logging.debug(sql.replace('?', '{}').format(min_cluster_prop))
    cursor.execute(sql, [min_cluster_prop])

    if exclude_clusters:
        clusters = (c for c in cursor if c[0] not in exclude_clusters)
    else:
        clusters = cursor

    grouped = itertools.groupby(clusters, operator.itemgetter(0))

    selected_clusters = set()
    futs = set()
    with futures.ThreadPoolExecutor(threads) as executor:
        for cluster_name, values in grouped:
            cluster_seq_names = sequences_hitting_cluster(
                deenurp_db, cluster_name)
            sample_weights = get_sample_weights(deenurp_db, cluster_seq_names)

            norm_sw = dict()
            for k, v in list(sample_weights.items()):
                norm_sw[k] = v / sample_total_weights[k]

            max_sample, max_weight = max(
                list(norm_sw.items()), key=operator.itemgetter(1))

            logging.info(
                'Cluster %s: Max hit by %s: %.3f%%, %d hits',
                cluster_name,
                max_sample,
                max_weight * 100,
                len(cluster_seq_names))

            cluster_refs = esl_sfetch_seqs(
                ref_fasta, cluster_members[cluster_name], ref_idx)

            # cluster_hit_seqs returns unicode: convert to string.
            query_seqs = esl_sfetch_seqs(
                fasta_file,
                (str(i) for i in cluster_seq_names),
                fa_idx)

            if max_weight < min_cluster_prop:
                msg = 'ID: {} max_weight {} < min_mass {}, skipping'
                msg = msg.format(cluster_name, max_weight, min_cluster_prop)
                logging.info(msg)
                continue

            selected_clusters.add(cluster_name)
            futs.add(executor.submit(
                select_sequences_for_cluster, cluster_refs, query_seqs,
                cluster_name=cluster_name,
                norm_sw=norm_sw,
                cluster_weight=sum(v[-1] for v in values),
                max_sample=max_sample,
                max_weight=max_weight,
                keep_leaves=refs_per_cluster))

        # Whitelist
        if include_clusters:
            for cluster in include_clusters:
                if cluster in selected_clusters:
                    msg = 'Sequences for whitelist cluster {} already selected'
                    logging.info(msg.format(cluster))
                    continue
                else:
                    cluster_refs = esl_sfetch_seqs(
                        ref_fasta, cluster_members[cluster], ref_idx)
                    futs.add(
                        executor.submit(
                            select_sequences_for_whitelist_cluster,
                            cluster_refs,
                            cluster,
                            keep_leaves=refs_per_cluster))

        while futs:
            try:
                done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
                futs = set(pending)
                for f in done:
                    if f.exception():
                        raise f.exception()

                    if exclude_sequences:
                        result = (r for r in f.result()
                                  if r.id not in exclude_sequences)
                    else:
                        result = f.result()

                    for ref in result:
                        yield ref
            except futures.TimeoutError:
                pass  # Keep waiting
            except:
                logging.exception("Caught error in child thread - exiting")
                executor.shutdown(False)
                raise
