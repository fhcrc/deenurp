#!/bin/bash

set -e

rm -rf output
mkdir output

BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}
$DEENURP filter-outliers \
	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
	 output/filtered.fasta --filtered-seqinfo output/filtered.seqinfo.csv
