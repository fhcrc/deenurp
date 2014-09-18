"""
Wrappers and context managers around external programs.
"""
import collections
import contextlib
import csv
import functools
import logging
import os
import os.path
import subprocess

from Bio import SeqIO
import peasel
from taxtastic.refpkg import Refpkg

from .util import as_fasta, ntf, tempdir, nothing, maybe_tempfile, which, require_executable, MissingDependencyError

DEFAULT_CMALIGN_THREADS = 1

"""Path to item in data directory"""
data_path = functools.partial(os.path.join, os.path.dirname(__file__), 'data')

"""16S bacterial covariance model"""
CM = data_path('RRNA_16S_BACTERIA.cm')


@contextlib.contextmanager
def as_refpkg(sequences, name='temp.refpkg', threads=None):
    """
    Context manager yielding a temporary reference package for a collection of aligned sequences.

    Builds a tree with FastTree, creates a reference package, yields.
    """
    sequences = list(sequences)
    with ntf(prefix='fasttree-', suffix='.log') as log_fp, \
         ntf(prefix='fasttree-', suffix='.tre') as tree_fp, \
         tempdir(prefix='refpkg') as refpkg_dir:

        log_fp.close()

        fasttree(sequences, log_path=log_fp.name, output_fp=tree_fp, gtr=True, threads=threads)
        tree_fp.close()

        rp = Refpkg(refpkg_dir(name), create=True)
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

def fasttree(sequences, output_fp, log_path=None, quiet=True, gtr=False,
        gamma=False, threads=None, prefix=None):

    executable = 'FastTreeMP' if threads and threads > 1 else 'FastTree'
    if executable == 'FastTreeMP' and not which('FastTreeMP'):
        executable = 'FastTree'
        logging.warn("Multithreaded FastTreeMP not found. Using FastTree")
    require_executable(executable)

    env = os.environ.copy()
    if threads:
        env['OMP_NUM_THREADS'] = str(threads)
    cmd = (prefix or []) + [executable, '-nt']
    for k, v in (('-gtr', gtr), ('-gamma', gamma), ('-quiet', quiet)):
        if v:
            cmd.append(k)
    if log_path is not None:
        cmd.extend(['-log', log_path])

    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=output_fp, stdin=subprocess.PIPE, env=env)
    count = SeqIO.write(sequences, p.stdin, 'fasta')
    assert count
    p.stdin.close()
    p.wait()
    if not p.returncode == 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)

def guppy_redup(placefile, redup_file, output):
    require_executable('guppy')
    cmd = ['guppy', 'redup', '-m', placefile, '-d', redup_file, '-o', output]
    logging.debug(' '.join(cmd))
    subprocess.check_call(cmd)

def pplacer(refpkg, alignment, posterior_prob=False, out_dir=None, threads=2, quiet=True):
    """
    Run pplacer on the provided refpkg
    """
    require_executable('pplacer')
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
        logging.debug(' '.join(cmd))
        subprocess.check_call(cmd, stdout=stdout)

    assert os.path.exists(jplace)

    return jplace

def rppr_min_adcl(jplace, leaves, algorithm='pam', posterior_prob=False, point_mass=True,
        always_include=None):
    """
    Run rppr min_adcl on the given jplace file, cutting to the given number of leaves
    Returns the names of the leaves *to remove*.
    """
    require_executable('rppr')
    cmd = ['rppr', 'min_adcl', '--algorithm', algorithm, jplace, '--leaves',
           str(leaves)]
    if point_mass:
        cmd.append('--point-mass')
    if posterior_prob:
        cmd.append('--pp')
    if always_include:
        cmd.extend(('--always-include', always_include))
    logging.debug(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

def rppr_min_adcl_tree(newick_file, leaves, algorithm='pam',
        always_include=None):
    """
    Run rppr min_adcl_tree on the given newick tree file, cutting to the given number of leaves.

    Returns the names of the leaves *to remove*.
    """
    cmd = ['rppr', 'min_adcl_tree', '--algorithm', algorithm, newick_file, '--leaves',
           str(leaves)]
    if always_include:
        cmd.extend(('--always-include', always_include))
    logging.debug(' '.join(cmd))
    output = subprocess.check_output(cmd)
    return output.splitlines()

def _require_cmalign_11(cmalign='cmalign'):
    """
    Check for cmalign version 1.1, raising an error if not found
    """
    version_str = 'INFERNAL 1.1'
    cmd = [cmalign, '-h']
    o = subprocess.check_output(cmd)
    if version_str not in o:
        msg = ('cmalign 1.1 not found. '
                'Expected {0} in output of "{1}", got:\n{2}').format(version_str, ' '.join(cmd), o)
        raise MissingDependencyError(msg)

def cmalign_files(input_file, output_file, cm=CM, cpu=DEFAULT_CMALIGN_THREADS):
    cmd = ['cmalign']
    require_executable(cmd[0])
    _require_cmalign_11(cmd[0])
    cmd.extend(['--noprob', '--dnaout'])
    if cpu is not None:
        cmd.extend(['--cpu', str(cpu)])
    cmd.extend(['-o', output_file, cm, input_file])
    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(p.stdout.read().strip())
    error = p.stderr.read().strip()
    if p.wait() != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(p.returncode, error)

def cmalign(sequences, output=None, cm=CM, cpu=DEFAULT_CMALIGN_THREADS):
    """
    Run cmalign
    """
    with as_fasta(sequences) as fasta, \
         maybe_tempfile(output, prefix='cmalign', suffix='.sto', dir='.') as tf:
        cmalign_files(fasta, tf.name, cm=cm, cpu=cpu)

        for sequence in SeqIO.parse(tf, 'stockholm'):
            yield sequence


def esl_sfetch(sequence_file, name_iter, output_fp, use_temp=False):
    """
    Fetch sequences named in name_iter from sequence_file, indexing if
    necessary, writing to output_fp.

    If ``use_temp`` is True, a temporary index is created and used.
    """

    if use_temp:
        with peasel.temp_ssi(sequence_file) as index:
            sequences = (index[i] for i in name_iter)
            count = peasel.write_fasta(sequences, output_fp)
    else:
        # if not os.path.exists(sequence_file + '.ssi'):
        #     logging.warning("No index exists for %s. creating.", sequence_file)
        try:
            peasel.create_ssi(sequence_file)
        except IOError:
            logging.warning("An index already exists for %s", sequence_file)

        index = peasel.open_ssi(sequence_file)
        sequences = (index[i] for i in name_iter)
        count = peasel.write_fasta(sequences, output_fp)

    return count


def load_tax_maps(fps, has_header=False):
    """
    Load tax maps from an iterable of file pointers
    """
    d = {}
    for fp in fps:
        reader = csv.reader(fp)
        if has_header:
            next(reader)  # Skip
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
    require_executable('dnaclust')
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
