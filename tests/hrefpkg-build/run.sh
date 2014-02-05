#!/bin/bash

set -e

BASE=../rdp_10_30_named1200bp_subset
rm -rf hrefpkg
mkdir hrefpkg
DEENURP=${DEENURP-../../deenurp.py}
$DEENURP hrefpkg-build --index-rank=family $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv --output-dir hrefpkg
