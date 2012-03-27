"""
Select reference sequences for inclusion
"""
import json
import os
import os.path
import itertools
import logging
import subprocess
import sqlite3
import tempfile

from romperroom import as_fasta, fasttree, as_tree_file
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

CM = data_path('bacteria16S_508_mod5.cm')

def cmalign(sequences, mpi_args=None):
    if not mpi_args:
        cmd = ['cmalign']
    else:
        cmd = ['mpirun'] + mpi_args + ['cmalign', '--mpi']
    cmd.extend(['--sub', '-1', '--dna', '--hbanded'])
    with as_fasta(sequences) as fasta, open(os.devnull) as devnull, \
         tempfile.NamedTemporaryFile(prefix='cmalign', suffix='.sto', dir='.') as tf:
        cmd.extend(('-o', tf.name))
        cmd.append(CM)
        cmd.append(fasta)
        subprocess.check_call(cmd, stdout=devnull)

        for sequence in SeqIO.parse(tf, 'stockholm'):
            yield sequence

def vorotree(tree, leaves, algorithm='full', query_seqs=None):
    cmd = ['rppr', 'vorotree', '--algorithm', algorithm, tree, '--leaves',
            str(leaves)]
    if query_seqs:
        cmd.extend(('--query-seqs', ','.join(query_seqs)))
    logging.info(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

def seqrecord(name, residues, **annotations):
    sr = SeqRecord(Seq(residues), name)
    sr.annotations.update(annotations)
    return sr

def sequence_extractor(rdp_con):
    """
    returns a function to extract sequences that match a sequence name, or are a
    type strain from one of the represented tax_ids
    """
    def _generate_in(coll):
        return '({0})'.format(','.join('?' for i in coll))

    rdp_con.row_factory = sqlite3.Row
    cursor = rdp_con.cursor()

    stmt = """SELECT * FROM sequences
WHERE name IN {0} OR (tax_id IN {1} AND is_type = 'type');"""

    def extract_seqs(sequence_ids, tax_ids):
        s = stmt.format(_generate_in(sequence_ids), _generate_in(tax_ids))
        logging.info(s)
        cursor.execute(s, list(sequence_ids) + list(tax_ids))
        for i in cursor:
            sr = seqrecord(i['name'], i['residues'], species_name=i['species_name'],
                lineage=json.loads(i['lineage']), tax_id=i['tax_id'], is_type=i['is_type'])
            yield sr

    return extract_seqs

def select_sequences_for_cluster(ref_seqs, query_seqs, keep_leaves=5):
    """
    Given a set of reference sequences and query sequences
    """
    if len(ref_seqs) <= keep_leaves:
        return [i.id for i in ref_seqs]

    c = itertools.chain(ref_seqs, query_seqs)
    query_taxa = [i.id for i in query_seqs]
    aligned = cmalign(c)
    tree = fasttree(aligned)
    assert len(tree.get_terminals())  == len(ref_seqs) + len(query_seqs)
    with as_tree_file(tree) as tf:
        prune_leaves = frozenset(vorotree(tf, keep_leaves, query_seqs=query_taxa))
    result = frozenset(i.id for i in ref_seqs) - prune_leaves
    assert len(result) == keep_leaves
    return result

def choose_references(deenurp_db, rdp_con, refs_per_cluster=5):
    extractor = sequence_extractor(rdp_con)
    total_weight = deenurp_db.total_weight()

    for cluster_id, seqs, hits in deenurp_db.hits_by_cluster():
        if not hits:
            logging.info("No hits for cluster %d", cluster_id)
            continue

        hit_names = [i['best_hit_name'] for i in hits]
        tax_ids = frozenset(i['tax_id'] for i in hits)
        seqs = [seqrecord(i['name'], i['residues'], weight=i['weight']) for i in seqs]
        ref_seqs = list(extractor(hit_names, tax_ids))
        keep = select_sequences_for_cluster(ref_seqs, seqs, refs_per_cluster)
        refs = [i for i in ref_seqs if i.id in keep]
        assert len(refs) == len(keep)
        for i in refs:
            yield i
