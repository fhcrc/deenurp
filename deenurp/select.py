"""
Select reference sequences for inclusion
"""
import itertools
import logging
import operator
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .wrap import cmalign, as_refpkg, as_fasta, tempdir, redupfile_of_seqs, \
                  voronoi, guppy_redup, pplacer, esl_sfetch

DEFAULT_THREADS = 12
CLUSTER_THRESHOLD = 0.999

def seqrecord(name, residues, **annotations):
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr

def _sequence_extractor(con):
    """
    returns a function to extract sequences that match a sequence name, or are a
    type strain from one of the represented tax_ids
    """
    cursor = con.cursor()
    refpath_cache = {}
    def reference_path(ref_id):
        try:
            return refpath_cache[ref_id]
        except KeyError:
            cursor.execute("""SELECT file_path FROM refs WHERE ref_id = ?""",
                    [ref_id])
            r = cursor.fetchone()[0]
            refpath_cache[ref_id] = r
            return r

    def extract_seqs(id_ref):
        """
        Extract sequences from id, ref_id pairs
        """
        key = operator.itemgetter(1)
        id_ref = sorted(id_ref, key=key)

        with tempfile.NamedTemporaryFile(prefix='refs', suffix='.fasta') as tf:
            for g, v in itertools.groupby(id_ref, key):
                names = [i for i, _ in v]
                ref_path = reference_path(g)
                esl_sfetch(ref_path, names, tf)

            tf.seek(0)
            for s in SeqIO.parse(tf, 'fasta'):
                yield s

    return extract_seqs

def _cluster(sequences, threshold=CLUSTER_THRESHOLD):
    sequences = sorted(sequences, key=lambda s: len(s), reverse=True)
    with as_fasta(sequences) as fp, tempfile.NamedTemporaryFile() as ntf:
        cmd = ['usearch', '-cluster', fp, '-seedsout', ntf.name, '-id',
                str(threshold), '-quiet', '-nowordcountreject']
        subprocess.check_call(cmd)
        r = list(SeqIO.parse(ntf, 'fasta'))
    logging.info("Clustered %d to %d", len(sequences), len(r))
    return r

def select_sequences_for_cluster(ref_seqs, query_seqs, keep_leaves=5,
        threads=DEFAULT_THREADS, mpi_args=None):
    """
    Given a set of reference sequences and query sequences
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
        prune_leaves = set(voronoi(placedir('redup.jplace'), keep_leaves))

    result = frozenset(i.id for i in ref_seqs) - prune_leaves

    assert len(result) == keep_leaves

    return result

def choose_references(deenurp_db, refs_per_cluster=5, candidates=30,
        threads=DEFAULT_THREADS, min_cluster_prop=0.0, mpi_args=None,
        cluster_factor=1):
    """
    Choose reference sequences from a search, choosing refs_per_cluster
    reference sequences for each nonoverlapping cluster.
    """
    extractor = _sequence_extractor(deenurp_db.con)
    total_weight = deenurp_db.total_weight()

    for cluster_id, seqs, hits in deenurp_db.hits_by_cluster(candidates,
            cluster_factor=cluster_factor):
        if not hits:
            logging.debug("No hits for cluster %d", cluster_id)
            continue
        cluster_weight = sum(seq['weight'] for seq in seqs)
        cluster_weight_prop = cluster_weight / total_weight

        logging.info("Cluster #%d: %d references; %d sequences "
                "[weight: %.2f/%.2f (%.1f%%)]", cluster_id, len(hits), len(seqs),
                cluster_weight, total_weight, cluster_weight_prop * 100)

        if cluster_weight_prop < min_cluster_prop:
            logging.info("Not enough mass in cluster #%d. Skipping.", cluster_id)
            continue

        # Increase the number of sequences to select for merged clusters by
        # cluster_factor
        cluster = deenurp_db.get_cluster(cluster_id)
        select_count = min(refs_per_cluster + cluster_factor * (cluster.count - 1),
                len(hits))
        if cluster.count > 1:
            logging.info("Selecting %d sequences for %d merged clusters.",
                    select_count, cluster.count)

        hit_names = frozenset((i['best_hit_name'], i['ref_id']) for i in hits)
        seqs = [seqrecord(i['name'], i['residues'], weight=i['weight']) for i in seqs]
        ref_seqs = list(extractor(hit_names))
        keep = select_sequences_for_cluster(ref_seqs, seqs, select_count,
                threads=threads, mpi_args=mpi_args)
        refs = [i for i in ref_seqs if i.id in keep]

        # Annotate
        for sequence in refs:
            sequence.annotations.update({'weight_prop': cluster_weight_prop,
                'cluster_id': cluster_id})

        assert len(refs) == len(keep)

        for i in refs:
            yield i
