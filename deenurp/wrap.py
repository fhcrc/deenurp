"""
Wrappers and context managers
"""
import collections
import contextlib
import csv
import logging
import os
import os.path
import subprocess

from Bio import SeqIO
from taxtastic.refpkg import Refpkg

from .util import as_fasta, ntf, tempdir, nothing, maybe_tempfile

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

CM = data_path('bacteria16S_508_mod5.cm')

@contextlib.contextmanager
def as_refpkg(sequences, threads=None):
    """
    Build a tree from sequences, generate a temporary reference package
    """
    sequences = list(sequences)
    with ntf(prefix='fast', suffix='.log') as log_fp, \
         ntf(prefix='fast', suffix='.tre') as tree_fp, \
         tempdir(prefix='refpkg') as refpkg_dir:

        log_fp.close()

        fasttree(sequences, log_fp.name, tree_fp, gtr=True, threads=threads)
        tree_fp.close()

        rp = Refpkg(refpkg_dir('temp.refpkg'))
        rp.update_metadata('locus', '')
        rp.update_phylo_model('FastTree', log_fp.name)
        rp.update_file('tree', tree_fp.name)

        # FASTA and Stockholm alignment
        with ntf(suffix='.fasta') as f:
            SeqIO.write(sequences, f, 'fasta')
            f.close()
            rp.update_file('aln_fasta', f.name)
        with ntf(suffix='.sto') as f:
            SeqIO.write(sequences, f, 'stockholm')
            f.close()
            rp.update_file('aln_sto', f.name)
        logging.debug("Reference package written to %s", rp.path)
        yield rp

@contextlib.contextmanager
def redupfile_of_seqs(sequences, **kwargs):
    with ntf(**kwargs) as tf:
        writer = csv.writer(tf, lineterminator='\n')
        rows = ((s.id, s.id, s.annotations.get('weight', 1.0)) for s in sequences)
        writer.writerows(rows)
        tf.flush()
        tf.close()
        yield tf.name

def fasttree(sequences, log_path, output_fp, quiet=True, gtr=False,
        gamma=False, threads=None, prefix=None):
    executable = 'FastTreeMP' if threads and threads > 1 else 'FastTree'
    env = os.environ.copy()
    if threads:
        env['OMP_NUM_THREADS'] = str(threads)
    cmd = (prefix or []) + [executable, '-nt', '-log', log_path]
    for k, v in (('-gtr', gtr), ('-gamma', gamma), ('-quiet', quiet)):
        if v:
            cmd.append(k)

    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=output_fp, stdin=subprocess.PIPE, env=env)
    count = SeqIO.write(sequences, p.stdin, 'fasta')
    assert count
    p.stdin.close()
    p.wait()
    if not p.returncode == 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)

def guppy_redup(placefile, redup_file, output):
    cmd = ['guppy', 'redup', '-m', placefile, '-d', redup_file, '-o', output]
    logging.info(' '.join(cmd))
    subprocess.check_call(cmd)

def pplacer(refpkg, alignment, posterior_prob=False, out_dir=None, threads=2, quiet=True):
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

    stdout = open(os.devnull, 'w') if quiet else nothing()

    with stdout:
        logging.info(' '.join(cmd))
        subprocess.check_call(cmd, stdout=stdout)

    assert os.path.exists(jplace)

    return jplace

def voronoi(jplace, leaves, algorithm='pam', posterior_prop=False, point_mass=True,
        always_include=None):
    """
    Run rppr voronoi on the given jplace file, cutting to the given number of leaves
    """
    cmd = ['rppr', 'voronoi', '--algorithm', algorithm, jplace, '--leaves',
           str(leaves)]
    if point_mass:
        cmd.append('--point-mass')
    if posterior_prop:
        cmd.append('--pp')
    if always_include:
        cmd.extend(('--always-include', always_include))
    logging.info(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

def cmalign_files(input_file, output_file, mpi_args=None, cm=CM,
        stdout=None):
    if mpi_args is None:
        cmd = ['cmalign']
    else:
        cmd = ['mpirun'] + mpi_args + ['cmalign', '--mpi']
    cmd.extend(['--sub', '-1', '--dna', '--hbanded'])
    cmd.extend(['-o', output_file, cm, input_file])
    logging.debug(' '.join(cmd))
    subprocess.check_call(cmd, stdout=stdout)

def cmalign(sequences, output=None, mpi_args=None, cm=CM):
    """
    Run cmalign

    If mpi_args is specified, run via mpirun
    """
    with as_fasta(sequences) as fasta, open(os.devnull) as devnull, \
         maybe_tempfile(output, prefix='cmalign', suffix='.sto', dir='.') as tf:
        cmalign_files(fasta, tf.name, mpi_args=mpi_args, stdout=devnull, cm=cm)

        for sequence in SeqIO.parse(tf, 'stockholm'):
            yield sequence

def esl_sfetch(sequence_file, name_iter, output_fp):
    """
    Fetch sequences named in name_iter from sequence_file, indexing if
    necessary, writing to output_fp.
    """
    if not os.path.exists(sequence_file + '.ssi'):
        logging.info("No index exists for %s. creating.", sequence_file)
        subprocess.check_call(['esl-sfetch', '--index', sequence_file])
    cmd = ['esl-sfetch', '-f', sequence_file, '-']
    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=output_fp)
    count = 0
    for name in name_iter:
        p.stdin.write('{0}\n'.format(name))
        count += 1
    p.stdin.close()
    p.wait()
    if not p.returncode == 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)

    return count

def load_tax_maps(fps, has_header=False):
    """
    Load tax maps from an iterable of file pointers
    """
    d = {}
    for fp in fps:
        reader = csv.reader(fp)
        if has_header:
            next(reader) # Skip
        for row in reader:
            name, taxid = row[:2]
            if name in d and taxid != d[name]:
                raise ValueError("Multiple tax_ids specified for {0}".format(name))
            d[name] = taxid
    return d


def dnaclust(fasta_file, similarity=0.99, centers=None, left_gaps_allowed=True,
         no_overlap=False, approximate=False):
    """
    Run DNA clust
    """
    Cluster = collections.namedtuple('Cluster', ['center', 'sequences'])
    cmd = ['dnaclust', '-s', str(similarity), '-i', fasta_file]
    if no_overlap:
        cmd.append('--no-overlap')
    if left_gaps_allowed:
        cmd.append('-l')
    if centers:
        cmd.extend(('-p', centers))
    if approximate:
        cmd.append('--approximate-filter')

    with ntf(prefix='dnaclust-') as tf:
        subprocess.check_call(cmd, stdout=tf)
        tf.seek(0)
        results = (i.strip().split() for i in tf)
        results = (Cluster(i[0], set(i)) for i in results)
        for i in results:
            yield i
