#!/bin/bash

set -e

rm -rf output

BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

out=output/cmalign
mkdir -p $out
time $DEENURP filter-outliers --aligner cmalign \
     $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
     $out/filtered.fasta --filtered-seqinfo $out/filtered.seqinfo.csv

out=output/muscle
mkdir -p $out
time $DEENURP filter-outliers --aligner muscle \
     $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
     $out/filtered.fasta --filtered-seqinfo $out/filtered.seqinfo.csv
