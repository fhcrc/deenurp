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

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import search, uclust
from concurrent import futures

from . import util, wrap
from .config import DEFAULT_THREADS
from .util import as_fasta, tempdir
from .wrap import cmalign, as_refpkg, redupfile_of_seqs, \
    rppr_min_adcl, guppy_redup, pplacer, esl_sfetch

CLUSTER_THRESHOLD = 0.998

def log_error(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            logging.exception("ERROR in %s", fn.__name__)
            raise
    return wrapper

def seqrecord(name, residues, **annotations):
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr

def _cluster(sequences, threshold=CLUSTER_THRESHOLD):
    """
    Cluster ``sequences`` at ``threshold``, returning seeds.
    """
    sequences = list(sequences)
    assert sequences
    with as_fasta(sequences) as fasta_name, tempfile.NamedTemporaryFile(prefix='uc-') as ntf:
        uclust.sort_and_cluster(fasta_name, ntf.name, pct_id=threshold,
                quiet=True, wordcountreject=False)
        ntf.seek(0)
        r = list(uclust.cluster_seeds(fasta_name, ntf))

    logging.debug("Clustered %d to %d", len(sequences), len(r))
    return r

@log_error
def select_sequences_for_cluster(ref_seqs, query_seqs, cluster_name,
                                 cluster_weight, max_sample, max_weight, norm_sw, keep_leaves=5):
    """
    Given a set of reference sequences and query sequences, select
    keep_leaves appropriate references.
    """
    logging.info('Cluster %s: Max sample abundance: %.3f%% of %s, %d hits',
                 cluster_name, max_weight * 100, max_sample, len(query_seqs))
    # Cluster
    ref_seqs = _cluster(ref_seqs)
    for ref in ref_seqs:
        ref.annotations.update({'cluster_name': cluster_name,
                                'max_weight': max_weight,
                                'mean_weight': sum(norm_sw.values()) / len(norm_sw)})
    if len(ref_seqs) <= keep_leaves:
        return ref_seqs

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
        prune_leaves = set(rppr_min_adcl(placedir('redup.jplace'), keep_leaves))

    result = frozenset(i.id for i in ref_seqs) - prune_leaves
    assert len(result) == keep_leaves

    refs = [i for i in ref_seqs if i.id in result]
    return refs

@log_error
def select_sequences_for_whitelist_cluster(ref_seqs, cluster_name, keep_leaves=5):
    """
    Selects a subset of ``ref_seqs`` using ``rppr min_adcl_tree``
    """
    logging.info("Whitelisted cluster %s", cluster_name)

    ref_seqs = _cluster(ref_seqs)
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
    return ref_ids - frozenset(prune)

def fetch_cluster_members(cluster_info_file, group_field):
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
        cursor.execute(sql, [sequence])
        res = cursor.fetchone()
        if not res:
            continue
        sample, weight = res
        weights[sample] += weight
    return dict(weights)

def sequences_hitting_cluster(con, cluster_name):
    sql = """SELECT DISTINCT sequences.name
        FROM sequences
        INNER JOIN best_hits USING (sequence_id)
        INNER JOIN ref_seqs USING(ref_id)
        WHERE cluster_name = ?
        ORDER BY sequences.name"""
    cursor = con.cursor()
    cursor.execute(sql, [cluster_name])
    return [name for name, in cursor]

def esl_sfetch_seqs(sequence_file, sequence_names, **kwargs):
    with tempfile.NamedTemporaryFile(prefix='esl', suffix='.fasta') as tf:
        esl_sfetch(sequence_file, sequence_names, tf, **kwargs)
        tf.seek(0)
        return list(SeqIO.parse(tf, 'fasta'))

def get_total_weight_per_sample(con):
    sql = """SELECT name, SUM(weight) AS weight
    FROM sequences_samples
      INNER JOIN samples USING (sample_id)
    GROUP BY name"""
    cursor = con.cursor()
    cursor.execute(sql)
    return dict(cursor)

def choose_references(deenurp_db, refs_per_cluster=5,
        threads=DEFAULT_THREADS, min_cluster_prop=0.0, whitelist=None):
    """
    Choose reference sequences from a search, choosing refs_per_cluster
    reference sequences for each nonoverlapping cluster.
    """
    whitelist = whitelist or set()
    params = search.load_params(deenurp_db)
    fasta_file = params['fasta_file']
    ref_fasta = params['ref_fasta']
    sample_total_weights = get_total_weight_per_sample(deenurp_db)
    cluster_members = fetch_cluster_members(params['ref_meta'], params['group_field'])

    # Iterate over clusters
    cursor = deenurp_db.cursor()

    # Select all clusters above cutoff
    cursor.execute("""
SELECT ref_seqs.cluster_name, sample_id, SUM(sequences_samples.weight) AS total_weight
FROM ref_seqs
    INNER JOIN best_hits USING (ref_id)
    INNER JOIN sequences USING (sequence_id)
    INNER JOIN sequences_samples USING (sequence_id)
GROUP BY ref_seqs.cluster_name, sample_id
HAVING SUM(sequences_samples.weight) > ?
ORDER BY ref_seqs.cluster_name ASC, SUM(sequences_samples.weight) DESC
""", [min_cluster_prop])
    grouped = itertools.groupby(cursor, operator.itemgetter(0))

    selected_clusters = set()
    futs = set()
    with futures.ThreadPoolExecutor(threads) as executor:
        for cluster_name, values in grouped:
            cluster_seq_names = sequences_hitting_cluster(deenurp_db, cluster_name)
            sample_weights = get_sample_weights(deenurp_db, cluster_seq_names)
            norm_sw = {k: v / sample_total_weights[k] for k, v in sample_weights.items()}
            max_sample, max_weight = max(norm_sw.items(), key=operator.itemgetter(1))

            logging.info('Cluster %s: Max hit by %s: %.3f%%, %d hits',
                    cluster_name, max_sample, max_weight * 100, len(cluster_seq_names))

            assert cluster_name in cluster_members
            cluster_refs = esl_sfetch_seqs(ref_fasta, cluster_members[cluster_name])

            # cluster_hit_seqs returns unicode: convert to string.
            query_seqs = esl_sfetch_seqs(fasta_file, (str(i) for i in cluster_seq_names),
                    use_temp=True)

            if max_weight < min_cluster_prop:
                logging.info("Skipping.")
                continue

            selected_clusters.add(cluster_name)
            futs.add(executor.submit(select_sequences_for_cluster,
                cluster_refs, query_seqs, cluster_name=cluster_name,
                norm_sw=norm_sw, cluster_weight=sum(v[-1] for v in values),
                max_sample=max_sample, max_weight=max_weight,
                keep_leaves=refs_per_cluster))
        # Whitelist
        for cluster in whitelist:
            if cluster in selected_clusters:
                logging.info("Sequences for whitelist cluster %s already selected", cluster)
                continue
            else:
                cluster_refs = esl_sfetch_seqs(ref_fasta, cluster_members[cluster])
                futs.add(executor.submit(select_sequences_for_whitelist_cluster,
                    cluster_refs, cluster, keep_leaves=refs_per_cluster))

        while futs:
            try:
                done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
                futs = set(pending)
                for f in done:
                    if f.exception():
                        raise f.exception()
                    for ref in f.result():
                        assert hasattr(ref, 'id')
                        yield ref
            except futures.TimeoutError:
                pass  # Keep waiting
            except:
                logging.exception("Caught error in child thread - exiting")
                executor.shutdown(False)
                raise
