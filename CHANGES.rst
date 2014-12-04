=====================
 Changes for deenurp
=====================

0.0.4
=====

 * ``filter-outliers`` can perform alignments using either cmalign or muscle


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