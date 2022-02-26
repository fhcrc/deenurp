import os
import sys
import subprocess

from setuptools import setup, find_packages, Command

datadir = 'deenurp/data'
version_file = f'{datadir}/version.txt'

subprocess.call(
    (f'mkdir -p {datadir} && '
     f'git describe --tags --dirty > {version_file}.tmp '
     f'&& mv {version_file}.tmp {version_file} '
     f'|| rm -f {version_file}.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

# import must follow 'git describe' command above to update version
from deenurp import __version__


class CheckVersion(Command):
    description = 'Confirm that the stored package version is correct'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        with open(version_file) as f:
            stored_version = f.read().strip()

        git_version = subprocess.check_output(
            ['git', 'describe', '--tags', '--dirty']).strip()

        assert stored_version == git_version
        print('the current version is', stored_version)


class run_audit(Command):

    """Audits source code using PyFlakes for following issues:
        - Names which are used but not defined or used before they are defined.
        - Names which are redefined without having been used.
    """
    description = "Audit source code with PyFlakes"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            import pyflakes.scripts.pyflakes as flakes
        except ImportError:
            print("Audit requires PyFlakes installed in your system.")
            sys.exit(-1)

        warns = 0
        # Define top-level directories
        dirs = 'deenurp',
        for dir in dirs:
            for root, _, files in os.walk(dir):
                for file in files:
                    if file != '__init__.py' and file.endswith('.py'):
                        warns += flakes.checkPath(os.path.join(root, file))
        if warns > 0:
            print("Audit finished with total %d warnings." % warns)
        else:
            print("No problems found in sourcecode.")


setup(
    name='deenurp',
    version=__version__,
    package_data={'deenurp': ['data/*', 'test/data/*']},
    entry_points={
        'console_scripts': {'deenurp = deenurp:main'}},
    cmdclass={'audit': run_audit, 'check_version': CheckVersion},
    test_suite='deenurp.test.suite',
    packages=find_packages(exclude=['tests']),
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'cython',
        'pandas',
        'scipy',
        'scikit-learn',
        'hdbscan',
        'biopython',
        'taxtastic',
        'seqmagick',
    ],
)
