=====================
 Changes for deenurp
=====================

0.3.2
=====

* Dependency updates and Python3 conversion bug fixes

0.3.1
=====

* building Docker image and test casese using Github Actions

0.3.0
=====

* Python 3 support and other dependency upgrades and bug fixes

0.2.7
=====

* removed peasel requirement

0.2.6
=====

* ``filter_outliers`` performance improvements

0.2.5
=========

* new ``deenurp filter_outliers --no-filter`` list of taxids argument [GH-65]
* fixed some Python Pandas warnings (GH: 50)
* fix error when query and reference sequences have name collision

0.2.4
=====

* update infernal to 1.1.2
* install cmsearch and cmcalibrate
* update taxtastic to 0.8.5
* update versions of various other dependencies

0.2.3
=====

* add seqmagick to requirements.txt
* create moinpoints in Docker image

0.2.2
=====

* fix version number in docker and singularity containers

0.2.1
=====

* bugfix: fix error parsing cmalign scores file
* add Singularity build instructions

0.2.0
=====

* add subcommand 'dereplicate_named'
* add new clustering method to filter_outliers (--cluster-type RobustSingleLinkage)
* update defaults for filter_outliers
* use vsearch 2.5.0
* fix error in fill_lonely raised by node with no parent
* added `deenurp orientate_sequences --out_notmatched_taxids`
* address edge case in select_references which too few refs
  representing a species remain after clustering at CLUSTER_THRESHOLD
* remove subcommand `tax2tree` (along with dependencies on t2t and cogent)
* remove subcommand `rdp_extract_genbank`
* add Dockerfile (thanks, Sam Minot!)
* update many dependencies

0.1.8
======

* ``deenurp ncbi_extract_genbank`` will die after 10 query retries

0.1.7
=====

* replace uclust with vsearch (uclust dropped as dependency)
* ``bootstrap.sh`` installs python requirements in specified order with ``--no-deps``
* installs vsearch 1.10.2

0.1.3
=====
* bumping pplacer support to v1.1.alpha17
* adding pandas support to versions 0.17.*
* ``deenurp orientate_sequences`` takes a training set (data/types.fasta) and aligns sequences reverse complementing when necessary (GH-34)
* travis caching is setup for pip installs (GH-22)
* fixed pep 440 versioning bug when installing and running unittests (GH-31)
* new ``deenurp select-references`` ``--include-clusters`` ``--exclude-clusters`` ``--exclude-sequences`` arguments
* new ``deenurp fill_lonely`` ``--include-taxid`` ``--exclude-taxid`` arguments

0.1.2
=====
* ``deenurp rdp_extract_sequences`` now updates old tax_ids instead of dropping them (GH-27)
* ``deenurp rdp_extract_sequences`` creates new column taxid_classified instead of dropping tax_ids of unclassified sequences (GH-28)
* function tax_of_genbank no longer returns None if string 'uncultured bacterium' is in the organism name (GH-29)
* expanding type strain designations to include more terms (GH-16)

0.1.1
=====

* ``gb2csv`` creates a csv of Genbank records and optionally the references (GH-23)
* using new container based Travis build (GH-24)

0.1.0
=====

* ``filter-outliers`` can calculate pairwise distances using cmalign, muscle, or vsearch
* ``filter-outliers`` can use hierarchical clustering for outlier detection (``--strategy cluster``)
* output of ``filter-outliers --detailed-seqinfo`` includes centroid, distance from centroid, etc
* add dependencies on pandas, scipy, and scikit-learn
* bootstrap script installs muscle and vsearch
* add subcommand ``deenurp pairwise-distances``
* version number defined using ``git describe --tags --dirty``

0.0.3
=====

* ``deenurp search-sequences`` learned option --search-threshold to
  allow more permissive database searches
* Suppress FastTree messages in the absence of an error

0.0.2
=====

* require Infernal 1.1
* update 16S alignment model (from https://github.com/rdpstaff/fungene_pipeline/blob/eab8ab3751da687b4d6dbd553f6a1d8261d98385/resources/RRNA_16S_BACTERIA/model.cm)
* add bin/boostrap.sh to create execution environment with installed dependencies
* more fine-grained control over threading/multiprocessing
