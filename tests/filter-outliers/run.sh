#!/bin/bash

set -e

rm -rf output

BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

for aligner in cmalign muscle vsearch; do
# for aligner in vsearch; do
    out=output/$aligner/radius
    mkdir -p $out
    time $DEENURP filter_outliers \
	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
	 $out/filtered.fasta --filtered-seqinfo $out/filtered.seqinfo.csv \
	 --detailed-seqinfo $out/filtered.details.csv \
	 --aligner $aligner \
	 --strategy radius

    out=output/$aligner/cluster
    mkdir -p $out
    time $DEENURP filter_outliers \
	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
	 $out/filtered.fasta --filtered-seqinfo $out/filtered.seqinfo.csv \
	 --detailed-seqinfo $out/filtered.details.csv \
	 --aligner $aligner \
	 --strategy cluster \
	 --distance-percentile 90

done
