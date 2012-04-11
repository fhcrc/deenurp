"""
Select reference sequences for inclusion
"""
import contextlib
import csv
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
CLUSTER_THRESHOLD = 0.998

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
        d = dict(id_ref)

        with tempfile.NamedTemporaryFile(prefix='refs', suffix='.fasta') as tf:
            for g, v in itertools.groupby(id_ref, key):
                names = [i for i, _ in v]
                ref_path = reference_path(g)
                esl_sfetch(ref_path, names, tf)

            tf.seek(0)
            for s in SeqIO.parse(tf, 'fasta'):
                s.annotations['ref'] = d[s.id]
                yield s

    return extract_seqs

def _cluster(sequences, threshold=CLUSTER_THRESHOLD):
    with as_fasta(sequences) as fp, tempfile.NamedTemporaryFile() as ntf:
        cmd = ['usearch', '-cluster', fp, '-seedsout', ntf.name, '-id',
                str(threshold),
                '-usersort',
                '-quiet', '-nowordcountreject']
        subprocess.check_call(cmd)
        r = frozenset(i.id for i in SeqIO.parse(ntf, 'fasta'))
    logging.info("Clustered %d to %d", len(sequences), len(r))
    return [i for i in sequences if i.id in r]

@contextlib.contextmanager
def _always_include(ref_seqs, keep_leaves):
    key = lambda x: x.annotations['ref']
    ref_seqs = sorted(ref_seqs, key=key)
    seqs = list(next(itertools.groupby(ref_seqs, key=key))[1])
    if len(seqs) < keep_leaves:
        with tempfile.NamedTemporaryFile(prefix='keep') as tf:
            for i in seqs:
                print >> tf, i.id
            tf.flush()
            yield ref_seqs, tf.name
    else:
        yield seqs, None

def select_sequences_for_cluster(ref_seqs, query_seqs, keep_leaves=5,
        threads=DEFAULT_THREADS, mpi_args=None):
    """
    Given a set of reference sequences and query sequences, select appropriate
    keep_leaves appropriate references.
    """
    # Cluster
    ref_seqs = _cluster(ref_seqs)
    if len(ref_seqs) <= keep_leaves:
        return [i.id for i in ref_seqs]

    with _always_include(ref_seqs, keep_leaves) as (ref_seqs, always):
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
            prune_leaves = set(voronoi(placedir('redup.jplace'), keep_leaves,
                              always_include=always))

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

        # Select sequences for the cluster, using `rppr voronoi`
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

def merge_meta(sequence_file, deenurp_db, output_fp):
    """
    Generate a merged metadata file for all sequences in sequence_file
    """
    seen_headers = set()
    headers = []
    sequence_ids = set(i.id for i in SeqIO.parse(sequence_file, 'fasta'))
    assert sequence_ids
    result = []

    cursor = deenurp_db.con.cursor()
    cursor.execute('SELECT DISTINCT meta_path FROM refs WHERE meta_path IS NOT NULL')
    meta_files = [i[0] for i in cursor]

    for meta_file in meta_files:
        with open(meta_file) as fp:
            reader = csv.DictReader(fp)
            h = set(reader.fieldnames)
            headers.extend(h - seen_headers)
            seen_headers |= set(reader.fieldnames)
            for i in reader:
                if i['seqname'] in sequence_ids:
                    sequence_ids.remove(i['seqname'])
                    result.append(i)

    # Add an empty line for everything else
    for i in sequence_ids:
        result.append({'seqname': i})

    writer = csv.DictWriter(output_fp, headers)
    writer.writeheader()
    writer.writerows(result)

    return len(result)
