"""
Utility functions
"""

import logging
import os
import os.path
import pkg_resources
import re
import subprocess


def version():
    """
    Depending on if git exists and at what version we can try these commands
    to get the git version from the git tags.  Second, if the package has been
    installed return the installed version.  Finally, if nothing else,
    return 0.0.0
    """

    version = None

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
        logging.debug(' '.join(cmd))
        git_re = 'v(?P<tag>[\d.]*)-?(?P<commit>[\d.]*)-?.*'
        try:
            git_ver = subprocess.check_output(cmd, stderr=devnull)
            git_search = re.search(git_re, git_ver)
            if git_search.group('commit') == '':
                version = git_search.group('tag')
            else:
                version = '{tag}.dev{commit}'.format(**git_search.groupdict())
            break
        except Exception as e:
            logging.debug(e)

    if version is None:
        try:
            """
            return version that was installed if available
            """
            version = pkg_resources.require("deenurp")[0].version
        except pkg_resources.DistributionNotFound as e:
            logging.debug(e)

    return version or '0.0.0'
