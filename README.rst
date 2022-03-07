=========
 DeeNuRP
=========

16S rRNA gene sequence curation and phylogenetic reference set creation

Installation
============

The Easy Way
------------

* confirm availability of necessary libraries to compile dependencies
  (on Ubuntu: ``sudo apt-get install gfortran libopenblas-dev liblapack-dev``)
* Install Python >= 3.8 or Python 3 Virtual Environment 
::

  % python3 -m venv bin-env
  % source bin-env/bin/activate
  % bin/bootstrap.sh


the `deenurp` executable should now be on your `$PATH`

The Hard Way
------------

See required system libraries above.

First, install binary dependencies:

* Python 3

    - pip, for installing python dependencies (http://www.pip-installer.org/)
    - Python packages:

        + Run ``pip install PACKAGE`` for every PACKAGE listed in requirements.txt, e.g.
          ``cat requirements.txt | xargs -n 1 pip install``

* vsearch (https://github.com/torognes/vsearch)
* Infernal version 1.1 (http://infernal.janelia.org/)
* pplacer suite (http://matsen.fhcrc.org/pplacer)
* FastTree 2 (http://www.microbesonline.org/fasttree/#Install)

Optional (for ``filter-outliers`` and ``pairwise-distances``):

* muscle (http://www.drive5.com/muscle/)

Finally, install::

    python setup.py install

The Docker Way
==============

Deenurp can be run from a Docker image which can be built locally from the Dockerfile
or pulled ``docker pull nghoffman/deenurp:v0.3.0``

De-novo reference set creation
==============================

Similarity-search based reference sequence selection

Running
-------

The ``deenurp`` package under the current directory provides to subcommands,
accessed via the script ``deenurp.py``, or the command ``deenurp`` if installed.

Subcommands fall into two general categories:

* Building a set of reference sequences for use in refpkg building
* Selecting sequences for a specific reference package

Creating a sequence set for refpkg building
-------------------------------------------

``deenurp filter-outliers``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Removes outlier sequences from a reference database

``deenurp expand-named``
~~~~~~~~~~~~~~~~~~~~~~~~

Expands poorly-represented names in a sequence file by similarity search

``deenurp cluster-refs``
~~~~~~~~~~~~~~~~~~~~~~~~

Cluster reference sequences, first by tax-id at a specified rank
(default: species), then by similarity for unnamed sequences or
sequences not classified to the desired rank.  Serves as input to
``search-sequences``.

Selecting sequences for a reference package
-------------------------------------------

``deenurp hrefpkg-build``
~~~~~~~~~~~~~~~~~~~~~~~~~

Builds a set of hierarchical reference packages.

``deenurp search-sequences``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Searches a set of sequences against a FASTA file containing possible
reference sequences.

This subcommand does searches sequences against a reference FASTA
file, saving the results and some metadata to a sqlite database for
use in ``select-references``

``deenurp select-references``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given the output of ``search-sequences``, ``select-references``
attempts to find a good set of reference sequences.

For each reference cluster with a minimal amount of sequences having
best hits to the cluster, (see ``cluster-refs``), selects a set number
of sequences to serve as references.

``deenurp fill-lonely``
~~~~~~~~~~~~~~~~~~~~~~~

Taxa who are the sole descendent of their parent can complicate
taxonomic classification.

The ``fill-lonely`` subcommand finds some company for these lonely
taxa.

``deenurp add-reps``
~~~~~~~~~~~~~~~~~~~~

Fetches sequences from a sequence file which match the taxtable for a
reference set at a given rank. Useful for adding type strains.

``deenurp tax2tree``
~~~~~~~~~~~~~~~~~~~~

Runs the ``tax2tree`` program on a reference package, updating the
``seq_info`` file.

Sequences whose lineage changes are relabeled. The prior ``tax_id`` is
added to the ``seq_info`` file in the reference package.

