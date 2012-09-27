from setuptools import setup, find_packages

install_requires = ['biopython>=1.58',
        'tax2tree==1.0',
        'cogent>=1.5.1',
        'taxtastic>=0.4.0']

setup(name='deenurp',
      version='0.0.1',
      package_data={'deenurp': ['data/*.cm']},
      entry_points={'console_scripts': {'deenurp = deenurp.scripts.deenurp:main'}},
      install_requires=install_requires,
      packages=find_packages())
