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
from packaging.version import Version
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
        tree_fp.close()
        fasttree(sequences, log_path=log_fp.name, output_fp=tree_fp.name,
                 gtr=True, threads=threads)

        rp = Refpkg(refpkg_dir(name), create=True)
        rp.update_metadata('locus', '')
        rp.update_phylo_model('FastTree', log_fp.name)
        rp.update_file('tree', tree_fp.name)

        # FASTA and Stockholm alignment
        with ntf('w', suffix='.fasta') as f:
            SeqIO.write(sequences, f, 'fasta')
            f.close()
            rp.update_file('aln_fasta', f.name)
        with ntf('w', suffix='.sto') as f:
            SeqIO.write(sequences, f, 'stockholm')
            f.close()
            rp.update_file('aln_sto', f.name)
        logging.debug("Reference package written to %s", rp.path)
        yield rp


@contextlib.contextmanager
def redupfile_of_seqs(sequences):
    with ntf('w') as tf:
        writer = csv.writer(tf, lineterminator='\n')
        rows = ((s.id, s.id, s.annotations.get('weight', 1.0)) for s in sequences)
        writer.writerows(rows)
        tf.flush()
        tf.close()
        yield tf.name


def fasttree(sequences, output_fp, log_path=None, quiet=True,
             gtr=False, gamma=False, threads=FASTTREE_THREADS, prefix=None):

    nseqs = len(sequences)
    if nseqs < 3:
        raise ValueError(
            f'at least 3 sequences are required but {nseqs} were provided')

    executable = 'FastTreeMP' if threads and threads > 1 else 'FastTree'
    if executable == 'FastTreeMP' and not which('FastTreeMP'):
        executable = 'FastTree'
        logging.warn("Multithreaded FastTreeMP not found. Using FastTree")
    require_executable(executable)

    env = os.environ.copy()
    if threads:
        env['OMP_NUM_THREADS'] = str(threads)

    with ntf('w', suffix='.fasta') as fasta:
        assert SeqIO.write(sequences, fasta, 'fasta')
        fasta.flush()

        cmd = (prefix or []) + [executable]
        opts = [('-gtr', gtr), ('-gamma', gamma), ('-quiet', quiet)]
        cmd.extend([k for k, v in opts if v])

        if log_path:
            cmd.extend(['-log', log_path])

        cmd.extend(['-out', output_fp, '-nt', fasta.name])
        logging.debug(' '.join(cmd))

        job = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if not job.returncode == 0:
            logging.error(job.stderr)
            raise subprocess.CalledProcessError(job.returncode, cmd)


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
    job = subprocess.run(cmd, capture_output=True, text=True)
    return job.stdout.strip().splitlines()


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
    o = subprocess.run([cmalign, '-h'], capture_output=True, text=True)
    if version_str not in o.stdout:
        msg = ('cmalign 1.1 not found. '
               'Expected {0} in output of "{1}", got:\n{2}').format(
                   version_str, ' '.join(cmd), o)
        raise MissingDependencyError(msg)


def cmalign_scores(text):
    """
    Parse stdout of cmalign into a data.frame
    """
    dtypes = {
        "idx": int,
        "seq_name": str,
        "length": int,
        "cm_from": int,
        "cm_to": int,
        "trunc": str,
        "bit_sc": float,
        "avg_pp": str,
        "band_calc": float,
        "alignment": float,
        "total": float,
        "mem": float
        }
    return pd.read_csv(
        StringIO(text),
        comment="#",
        dtype=dtypes,
        index_col='seq_name',
        names=dtypes.keys(),
        sep='\\s+'
        )


def cmalign_files(input_file, output_file, cm=CM, cpu=CMALIGN_THREADS):

    executable = 'cmalign'
    require_executable(executable)
    _require_cmalign_11(executable)
    cmd = executable + ' --noprob --dnaout'

    if cpu is not None:
        cmd += ' --cpu ' + str(cpu)
    cmd += ' -o {} {} {}'.format(output_file, cm, input_file)

    logging.debug(cmd)

    # TODO: preserve output files (input_file, output_file) on error
    # using shell=True to avoid Docker cmalign segmentation fault error (SIG 11)
    job = subprocess.run(
        cmd, capture_output=True, check=True, shell=True, text=True)

    output = job.stdout.strip()
    logging.debug(output)
    scores = cmalign_scores(output)

    return scores


def cmalign(sequences, output=None, cm=CM, cpu=CMALIGN_THREADS):
    """Run cmalign. If provided, saves output to the file path
    corresponding to file-like object `output`.

    """
    with as_fasta(sequences) as fasta, maybe_tempfile(
            output, mode='w+', prefix='cmalign', suffix='.sto', dir='.') as tf:

        cmalign_files(fasta, tf.name, cm=cm, cpu=cpu)
        for sequence in SeqIO.parse(tf, 'stockholm'):
            yield sequence


def _require_vsearch_version(vsearch=VSEARCH, version=VSEARCH_VERSION):
    """
    Check for vsearch with a version >= `version`
    """

    output = subprocess.run([vsearch, '--version'], capture_output=True, text=True)
    vsearch = re.search(r'^vsearch v(?P<vstr>\d+\.\d+\.[^_]+)', output.stderr)
    ver = vsearch.groupdict()['vstr']

    if Version(ver) < Version(version):
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
    job = subprocess.run(cmd, capture_output=True, text=True)
    logging.debug(job.stdout)

    if job.returncode != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(job.returncode, job.stderr)


def muscle_files(input_file, output_file, maxiters=MUSCLE_MAXITERS):
    cmd = [
        'muscle',
        '-in', input_file,
        '-out', output_file,
        # TODO: set value based on number of sequences?
        '-maxiters', str(maxiters),
    ]
    logging.debug(' '.join(cmd))
    require_executable(cmd[0])

    job = subprocess.run(cmd, capture_output=True, text=True)
    logging.debug(job.stdout)

    if job.returncode != 0:
        # TODO: preserve output files (input_file, output_file)
        raise subprocess.CalledProcessError(p.returncode, job.stderr)


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
    necessary, writing binary data to open file object output_fp.
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
