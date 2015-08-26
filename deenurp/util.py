"""
Utility functions
"""

import bz2
import contextlib
import functools
import gzip
import logging
import os
import os.path
import pandas
import pkg_resources
import re
import shutil
import subprocess
import sys
import time
import tempfile

from Bio import SeqIO


class Counter(object):

    """
    Count objects processed in iterable. By default, progress is written to
    stderr every 1000 items.
    """

    def __init__(self, iterable, stream=sys.stderr, report_every=0.3,
                 prefix=''):
        self._it = iter(iterable)
        self.count = 0
        self.stream = stream
        self.report_every = report_every
        self.prefix = prefix
        self.start = time.clock()
        self.last = 0

    def _report(self):
        if self.stream:
            msg = '{0}{1:15d} [{2:10.2f}s]\r'
            msg = msg.format(self.prefix, self.count, time.clock()-self.start)
            self.stream.write(msg)

    def __iter__(self):
        for i in self._it:
            yield i
            self.count += 1
            now = time.clock()
            if now - self.last > self.report_every:
                self._report()
                self.last = now


class SingletonDefaultDict(dict):

    """
    Dictionary-like object that returns the same value, regardless of key
    """

    def __init__(self, val=None):
        self.val = val

    def __getitem__(self, key):
        return self.val

    def __contains__(self, key):
        return True


def memoize(fn):
    cache = {}

    @functools.wraps(fn)
    def inner(*args):
        try:
            return cache[args]
        except KeyError:
            result = fn(*args)
            cache[args] = result
            return result
    inner.cache = cache
    return inner


def unique(iterable, key=lambda x: x):
    """
    Choose unique elements from iterable, using the value returned by `key` to
    determine uniqueness.
    """
    s = set()
    for i in iterable:
        k = key(i)
        if k not in s:
            s.add(k)
            yield i


@contextlib.contextmanager
def nothing(obj=None):
    """
    The least interesting context manager.
    """
    yield obj


@contextlib.contextmanager
def ntf(**kwargs):
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
def tempcopy(path, **kwargs):
    """
    Create a temporary copy of ``path``, available for the duration of the
    context manager
    """
    prefix, suffix = os.path.splitext(os.path.basename(path))
    a = {'prefix': prefix, 'suffix': suffix}
    a.update(kwargs)
    with open(path) as fp, ntf(**a) as tf:
        shutil.copyfileobj(fp, tf)
        tf.close()
        yield tf.name


@contextlib.contextmanager
def tempdir(**kwargs):
    """
    Create a temporary directory for the duration of the context manager,
    removing on exit.

    :returns: a partially applied os.path.join, with name of the temporary
    directory as the first argument

    Example:
    >>> with tempdir(prefix='rubbish-') as td:         # doctest: +SKIP
    ...     print "Directory is:", td()
    ...     print "Put some data in:", td('file1.txt')
    Directory is: /tmp/rubbish-5AQFpo
    Put some data in: /tmp/rubbish-5AQFpo/file1.txt
    """
    td = tempfile.mkdtemp(**kwargs)
    try:
        yield functools.partial(os.path.join, td)
    finally:
        shutil.rmtree(td)


@contextlib.contextmanager
def as_fasta(sequences, **kwargs):
    """
    Write sequences to a temporary FASTA file. returns the name
    """
    if 'suffix' not in kwargs:
        kwargs['suffix'] = '.fasta'
    with ntf(**kwargs) as tf:
        SeqIO.write(sequences, tf, 'fasta')
        tf.flush()
        tf.close()
        yield tf.name


@contextlib.contextmanager
def maybe_tempfile(obj=None, **kwargs):
    """
    Returns a tempfile for the duration of the contextmanager if obj is not
    provided, otherwise returns obj.
    """
    if obj is not None:
        yield obj
    else:
        with ntf(**kwargs) as tf:
            yield tf


@contextlib.contextmanager
def cd(path):
    """
    Change directory to `path` for the duration of the context manager
    """
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


def file_opener(mode='r'):
    """
    Returns a function that behaves similarly to ``open(...)``,
    but opens compressed files for certain matching extensions, currently
    ``.bz2`` is treated as bzip2-compression, and ``.gz`` is treated as gzip.
    """

    def open_file(f):
        out = None
        if f is sys.stdout or f is sys.stdin:
            out = f
        elif f == '-':
            out = sys.stdin if 'r' in mode else sys.stdout
        elif f.endswith('.bz2'):
            out = bz2.BZ2File(f, mode)
        elif f.endswith('.gz'):
            out = gzip.open(f, mode)
        else:
            out = open(f, mode)
        return out

    return open_file


def read_csv(fn, compression=None, **kwargs):
    """Read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.
    """

    if compression is not None or fn in [sys.stdin, '-']:
        compression = compression
    else:
        suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
        compression = compression or suffixes.get(os.path.splitext(fn)[-1])

    kwargs['compression'] = compression

    return pandas.read_csv(fn, **kwargs)


def which(executable_name, dirs=None):
    """
    Find an executable in dirs.
    If ``dirs`` is not specified, searches $PATH
    """
    if not dirs:
        dirs = os.environ['PATH'].split(os.pathsep)
    search_paths = (os.path.join(p, executable_name) for p in dirs)
    executable_paths = (i for i in search_paths
                        if os.path.exists(i) and os.access(i, os.EX_OK))
    try:
        return next(executable_paths)
    except StopIteration:
        return None


class MissingDependencyError(ValueError):
    pass


def require_executable(executable_name):
    if not which(executable_name):
        raise MissingDependencyError(executable_name)


def version():
    """
    Depending on if git exists and at what version we can try these commands
    to get the git version from the git tags.  Second, if the package has been
    installed return the installed version.  Finally, if nothing else,
    return 0.0.0
    """

    install_dir = os.path.dirname(__file__)
    git_cmds = (['git', '-C', install_dir, 'describe', '--tags'],  # >= 1.8.5
                ['git', 'describe', '--tags'])  # < 1.8.5
    devnull = open(os.devnull, 'w')

    """
    try the two git commands above ^
    git versions are organized as
    [tag]-[number of commits over tag]-[commit id]
    ex - v0.1.2-38-g6a8e5e3
    """
    for cmd in git_cmds:
        try:
            logging.debug(' '.join(cmd))
            git_re = 'v(?P<tag>[\d.]*)-?(?P<commit>[\d.]*)-?.*'
            git_ver = subprocess.check_output(cmd, stderr=devnull)
            git_search = re.search(git_re, git_ver)
            if git_search.group('commit') == '':
                return git_search.group('tag')
            else:
                return '{tag}.dev{commit} '.format(**git_search.groupdict())
        except Exception as e:
            logging.warn('{} {}'.format(type(e), e.message))

    try:
        """
        return version that was installed if available
        """
        return pkg_resources.require("csvpandas")[0].version
    except pkg_resources.DistributionNotFound as e:
        logging.warn(e)

    return '0.0.0'
