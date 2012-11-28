import functools
import os
import os.path

"""path to test file"""
data_path = functools.partial(os.path.join, os.path.dirname(__file__), 'data')

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

