"""
Wrappers and context managers
"""
import contextlib
import csv
import logging
import os
import os.path
import shutil
import subprocess
import tempfile

from Bio import SeqIO
from taxtastic.refpkg import Refpkg

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

CM = data_path('bacteria16S_508_mod5.cm')

@contextlib.contextmanager
def nothing(obj):
    """
    The least interesting context manager.
    """
    yield obj

@contextlib.contextmanager
def _ntf(**kwargs):
    """
    Near-clone of tempfile.NamedTemporaryFile, but the file is deleted when the
    context manager exits, rather than when it's closed.
    """
    kwargs['delete'] = False
    tf = tempfile.NamedTemporaryFile(**kwargs)
    try:
        with tf:
            yield tf
    finally:
        os.unlink(tf.name)

@contextlib.contextmanager
def tempdir(**kwargs):
    """
    Create a temporary directory for the duration of the context manager, removing
    on exit.
    """
    td = tempfile.mkdtemp(**kwargs)
    def p(*args):
        return os.path.join(td, *args)
    try:
        yield p
    finally:
        shutil.rmtree(td)

@contextlib.contextmanager
def as_fasta(sequences, **kwargs):
    """
    Write sequences to a temporary FASTA file. returns the name
    """
    if 'suffix' not in kwargs:
        kwargs['suffix'] = '.fasta'
    with _ntf(**kwargs) as tf:
        SeqIO.write(sequences, tf, 'fasta')
        tf.flush()
        tf.close()
        yield tf.name

@contextlib.contextmanager
def as_refpkg(sequences):
    """
    Build a tree from sequences, generate a temporary reference package
    """
    sequences = list(sequences)
    with _ntf(prefix='fast', suffix='.log') as log_fp, \
         _ntf(prefix='fast', suffix='.tre') as tree_fp, \
         tempdir(prefix='refpkg') as refpkg_dir:

        log_fp.close()

        fasttree(sequences, log_fp.name, tree_fp, gtr=True)
        tree_fp.close()

        rp = Refpkg(refpkg_dir('temp.refpkg'))
        rp.update_metadata('locus', '')
        rp.update_phylo_model('FastTree', log_fp.name)
        rp.update_file('tree', tree_fp.name)
        logging.info("Reference package written to %s", rp.path)
        yield rp

@contextlib.contextmanager
def redupfile_of_seqs(sequences, **kwargs):
    with _ntf(**kwargs) as tf:
        writer = csv.writer(tf, lineterminator='\n')
        rows = ((s.id, s.id, s.annotations.get('weight', 1.0)) for s in sequences)
        writer.writerows(rows)
        tf.flush()
        tf.close()
        yield tf.name

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

def guppy_redup(placefile, redup_file, output):
    cmd = ['guppy', 'redup', '-m', placefile, '-d', redup_file, '-o', output]
    logging.info(' '.join(cmd))
    subprocess.check_call(cmd)

def pplacer(refpkg, alignment, posterior_prob=True, out_dir=None, threads=2):
    """
    Run pplacer on the provided refpkg
    """
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

def voronoi(jplace, leaves, algorithm='full', posterior_prop=True, point_mass=True):
    """
    Run rppr voronoi on the given jplace file, cutting to the given number of leaves
    """
    cmd = ['rppr', 'voronoi', '--algorithm', algorithm, jplace, '--leaves',
           str(leaves)]
    if point_mass:
        cmd.append('--point-mass')
    if posterior_prop:
        cmd.append('--pp')
    logging.info(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

def cmalign(sequences, mpi_args=None):
    """
    Run cmalign

    If mpi_args is specified, run via mpirun
    """
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
