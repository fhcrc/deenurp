#!/bin/bash
BASE=../rdp_10_30_named1200bp_subset
mkdir -p hrefpkg
DEENURP=${DEENURP-../../deenurp.py}
$DEENURP hrefpkg-build --index-rank=family $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv --output-dir hrefpkg --threads 6
