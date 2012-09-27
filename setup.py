from setuptools import setup, find_packages, Command

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
        import os, sys
        try:
            import pyflakes.scripts.pyflakes as flakes
        except ImportError:
            print "Audit requires PyFlakes installed in your system."
            sys.exit(-1)

        warns = 0
        # Define top-level directories
        dirs = ('flask', 'examples', 'scripts')
        for dir in dirs:
            for root, _, files in os.walk(dir):
                for file in files:
                    if file != '__init__.py' and file.endswith('.py') :
                        warns += flakes.checkPath(os.path.join(root, file))
        if warns > 0:
            print "Audit finished with total %d warnings." % warns
        else:
            print "No problems found in sourcecode."

#install_requires = ['biopython>=1.58',
        #'tax2tree',
        #'cogent>=1.5.1',
        #'taxtastic>=0.4.0']

install_requires = []

setup(name='deenurp',
      version='0.0.1',
      package_data={'deenurp': ['data/*', 'test/data/*']},
      entry_points={'console_scripts': {'deenurp = deenurp.scripts.deenurp:main'}},
      install_requires=install_requires,
      cmdclass={'audit': run_audit},
      packages=find_packages())
