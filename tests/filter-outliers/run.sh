#!/bin/bash

set -e

rm -rf output

BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

for aligner in cmalign muscle vsearch; do
    out=output/$aligner
    mkdir -p $out
    time $DEENURP filter-outliers --aligner $aligner \
	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
	 $out/filtered.fasta --filtered-seqinfo $out/filtered.seqinfo.csv \
	 --detailed-seqinfo $out/filtered.details.csv
done
