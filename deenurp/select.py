"""
Select reference sequences for inclusion
"""
import csv
import contextlib
import json
import os
import os.path
import itertools
import logging
import subprocess
import shutil
import sqlite3
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from taxtastic.refpkg import Refpkg


def fasttree(sequences, log_path, output_fp, quiet=True, gtr=False, gamma=False):
    cmd = ['FastTree', '-nt', '-log', log_path]
    for k, v in (('-gtr', gtr), ('-gamma', gamma), ('-quiet', quiet)):
        if v:
            cmd.append(k)

    logging.info(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=output_fp, stdin=subprocess.PIPE)
    SeqIO.write(sequences, p.stdin, 'fasta')
    p.stdin.close()
    p.wait()
    if not p.returncode == 0:
        raise subprocess.CalledProcessError(p.returncode)

@contextlib.contextmanager
def tempdir(**kwargs):
    td = tempfile.mkdtemp(**kwargs)
    def p(*args):
        return os.path.join(td, *args)
    try:
        yield p
    finally:
        shutil.rmtree(td)

@contextlib.contextmanager
def as_fasta(sequences, **kwargs):
    if 'suffix' not in kwargs:
        kwargs['suffix'] = '.fasta'
    with tempfile.NamedTemporaryFile(**kwargs) as tf:
        SeqIO.write(sequences, tf, 'fasta')
        tf.flush()
        yield tf.name

@contextlib.contextmanager
def as_refpkg(sequences):
    """
    Build a tree from sequences, generate a temporary reference package
    """
    sequences = list(sequences)
    with tempfile.NamedTemporaryFile(prefix='fast', suffix='.log') as log_fp, \
         tempfile.NamedTemporaryFile(prefix='fast', suffix='.tre') as tree_fp, \
         tempdir(prefix='refpkg') as refpkg_dir:
        fasttree(sequences, log_fp.name, tree_fp, gtr=True)
        tree_fp.flush()

        rp = Refpkg(refpkg_dir('temp.refpkg'))
        rp.update_metadata('locus', '')
        rp.update_phylo_model('FastTree', log_fp.name)
        rp.update_file('tree', tree_fp.name)
        logging.info("Reference package written to %s", rp.path)
        yield rp

@contextlib.contextmanager
def redupfile_of_seqs(sequences, **kwargs):
    with tempfile.NamedTemporaryFile(**kwargs) as tf:
        writer = csv.writer(tf, lineterminator='\n')
        rows = ((s.id, s.id, s.annotations.get('weight', 1.0)) for s in sequences)
        writer.writerows(rows)
        tf.flush()
        yield tf.name

def guppy_redup(placefile, redup_file, output):
    cmd = ['guppy', 'redup', '-m', placefile, '-d', redup_file, '-o', output]
    logging.info(' '.join(cmd))
    subprocess.check_call(cmd)

def pplacer(refpkg, alignment, posterior_prob=True, out_dir=None, threads=2):
    cmd = ['pplacer', '-j', str(threads), '-c', refpkg, alignment]
    if posterior_prob:
        cmd.append('-p')
    if out_dir:
        cmd.extend(('--out-dir', out_dir))

    jplace = os.path.basename(os.path.splitext(alignment)[0]) + '.jplace'
    if out_dir:
        jplace = os.path.join(out_dir, jplace)

    logging.info(' '.join(cmd))
    subprocess.check_call(cmd)
    assert os.path.exists(jplace)
    return jplace

def voronoi(jplace, leaves, algorithm='full'):
    cmd = ['rppr', 'voronoi', '--algorithm', algorithm, jplace, '--leaves',
           str(leaves)]
    logging.info(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

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
    ref_ids = frozenset(i.id for i in ref_seqs)
    aligned = list(cmalign(c))
    with as_refpkg(i for i in aligned if i.id in ref_ids) as rp, \
         as_fasta(aligned) as fasta, \
         tempdir(prefix='jplace') as placedir, \
         redupfile_of_seqs(query_seqs) as redup_path:
        jplace = pplacer(rp.path, fasta, out_dir=placedir())
        # Redup
        guppy_redup(jplace, redup_path, placedir('redup.jplace'))
        prune_leaves = set(voronoi(placedir('redup.jplace'), keep_leaves))

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
