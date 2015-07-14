=====================
 Changes for deenurp
=====================

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
