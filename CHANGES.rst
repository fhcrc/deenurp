=====================
 Changes for deenurp
=====================

0.1.3
=========
* bumping pplacer support to v1.1.alpha17
* adding pandas support to versions 0.17.*
* ``deenurp orientate_sequences`` takes a training set (data/types.fasta) and aligns sequences reverse complementing when necessary (issue: 34)
* travis caching is setup for pip installs (issue: 22)
* fixed pep 440 versioning bug when installing and running unittests (issue: 31)

0.1.2
=====
* ``deenurp rdp_extract_sequences`` now updates old tax_ids instead of dropping them (issue: 27)
* ``deenurp rdp_extract_sequences`` creates new column taxid_classified instead of dropping tax_ids of unclassified sequences (issue: 28)
* function tax_of_genbank no longer returns None if string 'uncultured bacterium' is in the organism name (issue: 29)
* expanding type strain designations to include more terms (issue: 16)

0.1.1
=====

* ``gb2csv`` creates a csv of Genbank records and optionally the references (issue: 23)
* using new container based Travis build (issue: 24)

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
