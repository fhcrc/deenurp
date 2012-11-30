"""
Select reference sequences for inclusion
"""
import collections
import csv
import itertools
import logging
import operator
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import search, uclust
from .util import as_fasta, tempdir
from .wrap import cmalign, as_refpkg, redupfile_of_seqs, \
                  rppr_min_adcl, guppy_redup, pplacer, esl_sfetch

DEFAULT_THREADS = 12
CLUSTER_THRESHOLD = 0.998

def seqrecord(name, residues, **annotations):
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr

def _cluster(sequences, threshold=CLUSTER_THRESHOLD):
    """
    Cluster ``sequences`` at ``threshold``, returning seeds.
    """
    with as_fasta(sequences) as fasta_name, tempfile.NamedTemporaryFile(prefix='uc-') as ntf:
        uclust.sort_and_cluster(fasta_name, ntf.name, pct_id=threshold,
                quiet=True, wordcountreject=False)
        ntf.seek(0)
        r = list(uclust.cluster_seeds(fasta_name, ntf))

    logging.debug("Clustered %d to %d", len(sequences), len(r))
    return r

def select_sequences_for_cluster(ref_seqs, query_seqs, keep_leaves=5,
        threads=DEFAULT_THREADS, mpi_args=None):
    """
    Given a set of reference sequences and query sequences, select
    keep_leaves appropriate references.
    """
    # Cluster
    ref_seqs = _cluster(ref_seqs)
    if len(ref_seqs) <= keep_leaves:
        return [i.id for i in ref_seqs]

    c = itertools.chain(ref_seqs, query_seqs)
    ref_ids = frozenset(i.id for i in ref_seqs)
    aligned = list(cmalign(c, mpi_args=mpi_args))
    with as_refpkg((i for i in aligned if i.id in ref_ids), threads=threads) as rp, \
             as_fasta(aligned) as fasta, \
             tempdir(prefix='jplace') as placedir, \
             redupfile_of_seqs(query_seqs) as redup_path:

        jplace = pplacer(rp.path, fasta, out_dir=placedir(), threads=threads)
        # Redup
        guppy_redup(jplace, redup_path, placedir('redup.jplace'))
        prune_leaves = set(rppr_min_adcl(placedir('redup.jplace'), keep_leaves))

    result = frozenset(i.id for i in ref_seqs) - prune_leaves

    assert len(result) == keep_leaves

    return result

def fetch_cluster_members(cluster_info_file, group_field):
    d = collections.defaultdict(list)
    with open(cluster_info_file) as fp:
        r = csv.DictReader(fp)
        for i in r:
            d[i[group_field]].append(i['seqname'])
    return d

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
    return weights

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
        threads=DEFAULT_THREADS, min_cluster_prop=0.0, mpi_args=None):
    """
    Choose reference sequences from a search, choosing refs_per_cluster
    reference sequences for each nonoverlapping cluster.
    """
    params = search.load_params(deenurp_db)
    fasta_file = params['fasta_file']
    ref_fasta = params['ref_fasta']
    sample_total_weights = get_total_weight_per_sample(deenurp_db)
    cluster_members = fetch_cluster_members(params['ref_meta'], params['group_field'])

    # Iterate over clusters
    cursor = deenurp_db.cursor()
    cursor.execute("""SELECT cluster_name, total_weight
            FROM vw_cluster_weights
            ORDER BY total_weight DESC""")

    for cluster_name, cluster_weight in cursor:
        cluster_seq_names = sequences_hitting_cluster(deenurp_db, cluster_name)
        sample_weights = get_sample_weights(deenurp_db, cluster_seq_names)
        norm_sw = {k: v / sample_total_weights[k] for k, v in sample_weights.items()}
        max_sample, max_weight = max(norm_sw.items(), key=operator.itemgetter(1))

        logging.info('Cluster %s: Max hit by %s: %.3f%%, %d hits',
                cluster_name, max_sample, max_weight * 100, len(cluster_seq_names))

        cluster_refs = esl_sfetch_seqs(ref_fasta, cluster_members[cluster_name])
        # sequences_hitting_cluster returns unicode: convert to string.
        query_seqs = esl_sfetch_seqs(fasta_file, (str(i) for i in cluster_seq_names),
                use_temp=True)

        if max_weight < min_cluster_prop:
            logging.info("Skipping.")
            continue

        ref_names = select_sequences_for_cluster(cluster_refs, query_seqs, mpi_args=mpi_args,
                keep_leaves=refs_per_cluster, threads=threads)
        refs = (i for i in cluster_refs if i.id in ref_names)
        for ref in refs:
            ref.annotations.update({'cluster_name': cluster_name,
                'max_weight': max_weight,
                'mean_weight': sum(norm_sw.values())/len(norm_sw)})
            yield ref
