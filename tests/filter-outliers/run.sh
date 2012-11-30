#!/bin/bash
BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}
$DEENURP filter-outliers $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv filtered.fasta --filtered-seqinfo filtered.seqinfo.csv
