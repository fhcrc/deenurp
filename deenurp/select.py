"""
Select reference sequences for inclusion
"""
import collections
import csv
import functools
import itertools
import logging
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import search, uclust
from concurrent import futures

from .util import as_fasta, tempdir
from .wrap import cmalign, as_refpkg, redupfile_of_seqs, \
                  rppr_min_adcl, guppy_redup, pplacer, esl_sfetch

DEFAULT_THREADS = 12
CLUSTER_THRESHOLD = 0.998

def log_error(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            logging.exception("ERROR in %s", fn.name)
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
    with as_fasta(sequences) as fasta_name, tempfile.NamedTemporaryFile(prefix='uc-') as ntf:
        uclust.sort_and_cluster(fasta_name, ntf.name, pct_id=threshold,
                quiet=True, wordcountreject=False)
        ntf.seek(0)
        r = list(uclust.cluster_seeds(fasta_name, ntf))

    logging.debug("Clustered %d to %d", len(sequences), len(r))
    return r

@log_error
def select_sequences_for_cluster(ref_seqs, query_seqs, cluster_name,
        cluster_weight, total_weight, keep_leaves=5):
    """
    Given a set of reference sequences and query sequences, select
    keep_leaves appropriate references.
    """
    logging.info('Cluster %s: %.3f%%, %d hits', cluster_name,
            cluster_weight / total_weight * 100, len(query_seqs))
    # Cluster
    ref_seqs = _cluster(ref_seqs)
    if len(ref_seqs) <= keep_leaves:
        return ref_seqs

    c = itertools.chain(ref_seqs, query_seqs)
    ref_ids = frozenset(i.id for i in ref_seqs)
    aligned = list(cmalign(c, mpi_args=None))
    with as_refpkg((i for i in aligned if i.id in ref_ids), threads=1) as rp, \
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
    for ref in refs:
        ref.annotations.update({'cluster_name': cluster_name,
            'weight_prop': cluster_weight/total_weight})
    return refs

def fetch_cluster_members(cluster_info_file, group_field):
    d = collections.defaultdict(list)
    with open(cluster_info_file) as fp:
        r = csv.DictReader(fp)
        for i in r:
            d[i[group_field]].append(i['seqname'])
    return d

def cluster_hit_seqs(con, cluster_name):
    sql = '''SELECT DISTINCT sequences.name, weight
        FROM sequences
        INNER JOIN best_hits USING (sequence_id)
        INNER JOIN ref_seqs USING(ref_id)
        WHERE cluster_name = ?'''
    cursor = con.cursor()
    cursor.execute(sql, [cluster_name])
    return list(cursor)

def esl_sfetch_seqs(sequence_file, sequence_names, **kwargs):
    with tempfile.NamedTemporaryFile(prefix='esl', suffix='.fasta') as tf:
        esl_sfetch(sequence_file, sequence_names, tf, **kwargs)
        tf.seek(0)
        return list(SeqIO.parse(tf, 'fasta'))

def get_total_weight(con):
    sql = 'SELECT SUM(weight) FROM sequences'
    cursor = con.cursor()
    cursor.execute(sql)
    return cursor.fetchone()[0]

def choose_references(deenurp_db, refs_per_cluster=5,
        threads=DEFAULT_THREADS, min_cluster_prop=0.0):
    """
    Choose reference sequences from a search, choosing refs_per_cluster
    reference sequences for each nonoverlapping cluster.
    """
    params = search.load_params(deenurp_db)
    fasta_file = params['fasta_file']
    ref_fasta = params['ref_fasta']
    total_weight = get_total_weight(deenurp_db)
    cluster_members = fetch_cluster_members(params['ref_meta'], params['group_field'])

    # Iterate over clusters
    cursor = deenurp_db.cursor()
    cursor.execute('''SELECT cluster_name, total_weight
            FROM vw_cluster_weights
            ORDER BY total_weight DESC''')

    futs = set()
    with futures.ThreadPoolExecutor(threads) as executor:
        for cluster_name, cluster_weight in cursor:
            cluster_seq_names = dict(cluster_hit_seqs(deenurp_db, cluster_name))
            cluster_refs = esl_sfetch_seqs(ref_fasta, cluster_members[cluster_name])
            weight_prop = cluster_weight / total_weight

            # Annotate with cluster information
            for ref in cluster_refs:
                ref.annotations.update({'cluster_name': cluster_name,
                                        'weight_prop': weight_prop })

            # cluster_hit_seqs returns unicode: convert to string.
            query_seqs = esl_sfetch_seqs(fasta_file, (str(i) for i in cluster_seq_names),
                    use_temp=True)
            for i in query_seqs:
                i.annotations['weight'] = cluster_seq_names[i.id]

            if cluster_weight / total_weight < min_cluster_prop:
                logging.info("Skipping cluster %s. Total weight: %.3f%%",
                        cluster_name, cluster_weight / total_weight * 100)
                break

            futs.add(executor.submit(select_sequences_for_cluster,
                cluster_refs, query_seqs,
                cluster_name=cluster_name, cluster_weight=cluster_weight,
                total_weight=total_weight, keep_leaves=refs_per_cluster))

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
                pass # Keep waiting
            except:
                logging.exception("Caught error in child thread - exiting")
                executor.shutdown(False)
                raise
