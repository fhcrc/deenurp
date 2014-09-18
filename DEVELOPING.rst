====================
 Developing deenurp
====================

creating a release
==================

- run tests!
- update the version number in ``setup.py``
- update CHANGES.txt
- git commit
- create a git tag for the release, using format 'vmajor.minor.release' (eg, v0.2.0)
- ``git push origin master``
- ``git push --tags``
- update PyPi (see below)


PyPi
====

Make sure to have updated the version string in ``setup.py`` first.

If you have not done so create a ~/.pypirc file::

  python setup.py register

Proceed to build and upload::

  python setup.py clean
  python setup.py sdist bdist_wheel
  twine upload dist/*
