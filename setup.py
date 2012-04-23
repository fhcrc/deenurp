from setuptools import setup, find_packages

setup(name='deenurp',
      version='0.0.1',
      package_data={'deenurp': ['data/*.cm']},
      entry_points={'console_scripts': {'deenurp = deenurp.scripts.deenurp:main'}},
      packages=find_packages())
