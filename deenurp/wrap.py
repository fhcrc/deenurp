"""
Wrappers and context managers around external programs.
"""
import contextlib
import csv
import functools
import logging
import os
import os.path
import subprocess
import re
from distutils.version import LooseVersion
from io import StringIO

import pandas as pd

from Bio import SeqIO
from taxtastic.refpkg import Refpkg

from .util import (as_fasta, ntf, tempdir, nothing, maybe_tempfile,
                   which, require_executable, MissingDependencyError)

CMALIGN_THREADS = 4

MUSCLE_MAXITERS = 2

VSEARCH = 'vsearch'
VSEARCH_VERSION = '2.0.3'
VSEARCH_IDDEF = 2
VSEARCH_THREADS = 2

FASTTREE_THREADS = 4

"""Path to item in data directory"""
data_path = functools.partial(os.path.join, os.path.dirname(__file__), 'data')

"""16S bacterial covariance model"""
CM = data_path('RRNA_16S_BACTERIA.cm')


@contextlib.contextmanager
def as_refpkg(sequences, name='temp.refpkg', threads=FASTTREE_THREADS):
    """Context manager yielding a temporary reference package for a
    collection of aligned sequences.

    Builds a tree with FastTree, creates a reference package, yields.

    """
    sequences = list(sequences)
    with ntf(prefix='fasttree-', suffix='.log') as log_fp, \
         ntf(prefix='fasttree-', suffix='.tre') as tree_fp, \
         tempdir(prefix='refpkg') as refpkg_dir:

        log_fp.close()

        fasttree(sequences, log_path=log_fp.name, output_fp=tree_fp, gtr=True,
                 threads=threads)
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


def fasttree(sequences, output_fp, log_path=None, quiet=True,
             gtr=False, gamma=False, threads=FASTTREE_THREADS, prefix=None):

    if len(sequences) < 3:
        raise ValueError(
            'at least 3 sequences are required but {} were provided'.format(len(sequences)))

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

    with ntf() as stderr:
        p = subprocess.Popen(cmd, stdout=output_fp, stdin=subprocess.PIPE,
                             stderr=stderr, env=env)

        count = SeqIO.write(sequences, p.stdin, 'fasta')
        assert count
        p.stdin.close()
        p.wait()
        if not p.returncode == 0:
            stderr.seek(0)
            logging.error(stderr.read())
            raise subprocess.CalledProcessError(p.returncode, cmd)


def guppy_redup(placefile, redup_file, output):
    require_executable('guppy')
    cmd = ['guppy', 'redup', '-m', placefile, '-d', redup_file, '-o', output]
    logging.debug(' '.join(cmd))
    subprocess.check_call(cmd)


def pplacer(refpkg, alignment, posterior_prob=False, out_dir=None,
            threads=2, quiet=True):
    """Run pplacer on the provided refpkg

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


def rppr_min_adcl(jplace, leaves, algorithm='pam', posterior_prob=False,
                  point_mass=True, always_include=None):
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


def rppr_min_adcl_tree(newick_file, leaves, algorithm='pam', always_include=None):
    """Run rppr min_adcl_tree on the given newick tree file, cutting to
    the given number of leaves.

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
               'Expected {0} in output of "{1}", got:\n{2}').format(
                   version_str, ' '.join(cmd), o)
        raise MissingDependencyError(msg)


def cmalign_scores(text):
    """
    Parse stdout of cmalign into a data.frame
    """

    header_rexp = re.compile(r'^#\s+idx')

    lines = []
    for line in text.splitlines():
        if header_rexp.search(line):
            line = ' ' + line[1:].replace(' (Mb)', '')
            # replace single spaces
            line = re.sub(r'(?<! ) (?! )', '_', line)
        elif line.startswith('#'):
            continue
        lines.append(line)

    buf = StringIO('\n'.join(lines))
    tab = pd.read_fwf(buf, index_col=1)
    return tab


def cmalign_files(input_file, output_file, cm=CM, cpu=CMALIGN_THREADS):
    cmd = ['cmalign']
    require_executable(cmd[0])
    _require_cmalign_11(cmd[0])
    cmd.extend(['--noprob', '--dnaout'])
    if cpu is not None:
        cmd.extend(['--cpu', str(cpu)])
    cmd.extend(['-o', output_file, cm, input_file])
    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.stdout.read().strip()
    logging.debug(output)

    scores = cmalign_scores(output)

    error = p.stderr.read().strip()
    if p.wait() != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(p.returncode, error)

    return scores


def cmalign(sequences, output=None, cm=CM, cpu=CMALIGN_THREADS):
    """
    Run cmalign
    """
    with as_fasta(sequences) as fasta, maybe_tempfile(
            output, prefix='cmalign', suffix='.sto', dir='.') as tf:

        cmalign_files(fasta, tf.name, cm=cm, cpu=cpu)

        for sequence in SeqIO.parse(tf, 'stockholm'):
            yield sequence


def _require_vsearch_version(vsearch=VSEARCH, version=VSEARCH_VERSION):
    """
    Check for vsearch with a version >= `version`
    """

    cmd = [vsearch, '--version']
    p = subprocess.Popen(
        cmd,
        stderr=subprocess.PIPE,
        stdout=open(os.devnull, 'w'))
    __, stderr = p.communicate()
    vsearch = re.search(r'^vsearch v(?P<vstr>\d+\.\d+\.[^_]+)', stderr)
    ver = vsearch.groupdict()['vstr']

    if LooseVersion(ver) < LooseVersion(version):
        raise MissingDependencyError(
            'vsearch version >= v{} is required, got v{}'.format(version, ver))


def vsearch_allpairs_files(input_file, output_file, executable=VSEARCH,
                           threads=VSEARCH_THREADS, iddef=VSEARCH_IDDEF):
    """Use vsearch to calculate all pairwise distances.

    """

    require_executable(executable)
    _require_vsearch_version()

    cmd = [executable,
           '--allpairs_global', input_file,
           '--strand', 'plus',
           '--qmask', 'none',
           '--id', '0',
           '--threads', str(threads),
           '--iddef', str(iddef),
           '--blast6out', output_file]

    logging.info(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(p.stdout.read().strip())
    error = p.stderr.read().strip()
    if p.wait() != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(p.returncode, error)


def muscle_files(input_file, output_file, maxiters=MUSCLE_MAXITERS):
    cmd = ['muscle']
    require_executable(cmd[0])

    cmd.extend(['-in', input_file])
    cmd.extend(['-out', output_file])

    # TODO: set value based on number of sequences?
    cmd.extend(['-maxiters', str(maxiters)])

    logging.debug(' '.join(cmd))
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(p.stdout.read().strip())
    error = p.stderr.read().strip()
    if p.wait() != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(p.returncode, error)


def read_seq_file(sequence_file):
    """
    Reads a fasta file and records the binary offsets of each sequence
    """
    fa_idx = {}
    with open(sequence_file, 'rb') as sf:
        offset = 0
        name = None
        for line in sf:
            line_len = len(line)
            if line.startswith(b'>'):
                if name is not None:
                    fa_idx[name].append(offset)
                name = line[1:].split()[0].strip().decode()
                fa_idx[name] = [offset]
            offset += line_len
        # last sequence
        if name:
            fa_idx[name].append(offset)
    return fa_idx


def esl_sfetch(sequence_file, name_iter, output_fp, fa_idx):
    """
    Fetch sequences named in name_iter from sequence_file, indexing if
    necessary, writing to output_fp.
    """
    count = 0
    with open(sequence_file, 'rb') as fi:
        for name in name_iter:
            indices = fa_idx[name]
            fi.seek(indices[0])
            output_fp.write(fi.read(indices[1] - indices[0]))
            count += 1
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
