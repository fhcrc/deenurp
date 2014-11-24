====================
 Developing deenurp
====================

running tests
=============

Unit tests can be run like this::

  python setup.py test

and functional tests like this::

  tests/run.sh


creating a release
==================

We're using a feature branch workfow whenever feasible. So when a
feature is complete:

- run tests!
- update the version number in ``setup.py``
- update CHANGES.txt
- make a final commit
- merge into master and run tests again
- create a git tag for the release, using format 'vmajor.minor.release' (eg, v0.2.0)
- ``git push origin master``
- ``git push --tags``
- update PyPi (see below)


PyPi
----

(not on PyPi yet...)


Make sure to have updated the version string in ``setup.py`` first.

If you have not done so create a ~/.pypirc file::

  python setup.py register

Proceed to build and upload::

  python setup.py clean
  python setup.py sdist bdist_wheel
  twine upload dist/*
