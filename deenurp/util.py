"""
Utility functions
"""

import bz2
import contextlib
import functools
import gzip
import itertools
import os
import os.path
import shutil
import sys
import time
import tempfile

from Bio import SeqIO


def apply_df_status(func, df, msg=''):
    """

    """
    tmp_column = 'index_number'
    row_count = float(len(df))
    df[tmp_column] = range(int(row_count))
    msg += ' {:.0%}\r'

    def apply_func(item, msg):
        sys.stderr.write(msg.format(item[tmp_column] / row_count))
        return func(item)

    df = df.apply(apply_func, args=[msg], axis=1)
    return df.drop(tmp_column, axis=1)


class Counter(object):

    """
    Count objects processed in iterable. By default, progress is written to
    stderr every 0.3 seconds of work.
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
def ntf(*args, **kwargs):
    """
    Near-clone of tempfile.NamedTemporaryFile, but the file is deleted when the
    context manager exits, rather than when it's closed.
    """
    kwargs['delete'] = False
    tf = tempfile.NamedTemporaryFile(*args, **kwargs)
    try:
        with tf:
            yield tf
    finally:
        os.unlink(tf.name)


@contextlib.contextmanager
def tempcopy(path):
    """Create a temporary copy of ``path``, available for the duration of
    the context manager

    """

    prefix, suffix = os.path.splitext(os.path.basename(path))
    with ntf(prefix=prefix, suffix=suffix) as tf:
        tf.close()
        shutil.copyfile(path, tf.name)
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
    with ntf('w+', **kwargs) as tf:
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


def file_opener(mode='r', buffering=-1):
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
            out = bz2.BZ2File(f, mode=mode, buffering=buffering)
        elif f.endswith('.gz'):
            out = gzip.open(f, mode=mode)
        else:
            out = open(f, mode=mode, buffering=buffering)
        return out

    return open_file


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


def chunker(iterable, n, fillvalue=None):
    """
    Continuously chunk an iterator n items at a time
    """

    while True:
        chunk = list(itertools.islice(iterable, n))
        if chunk:
            yield chunk
        else:
            return
