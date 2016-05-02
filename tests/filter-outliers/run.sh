#!/bin/bash

set -e

rm -rf output
rm -f seqs.fasta.ssi

BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

# for aligner in vsearch; do
for aligner in cmalign muscle vsearch; do
    echo $aligner

    out=output/$aligner/radius
    mkdir -p $out

    time $DEENURP filter_outliers \
    	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
    	 $out/filtered.fasta \
    	 --filtered-seqinfo $out/filtered.seqinfo.csv \
    	 --detailed-seqinfo $out/filtered.details.csv \
    	 --aligner $aligner \
    	 --strategy radius

    out=output/$aligner/cluster
    mkdir -p $out

    time $DEENURP filter_outliers \
    	 $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
    	 $out/filtered.fasta \
    	 --filtered-seqinfo $out/filtered.seqinfo.csv \
    	 --detailed-seqinfo $out/filtered.details.csv \
    	 --aligner $aligner \
    	 --strategy cluster \
    	 --distance-percentile 90

    time $DEENURP filter_outliers \
      $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv \
      $out/filtered2.fasta \
      --previous-details $out/filtered.details.csv \
      --filtered-seqinfo $out/filtered2.seqinfo.csv \
      --detailed-seqinfo $out/filtered2.details.csv \
      --aligner $aligner \
      --strategy cluster \
      --distance-percentile 90

    diff -q \
      <(grep '>' $out/filtered.fasta | sort) \
      <(grep '>' $out/filtered2.fasta | sort)

    diff -q \
      <(grep '>' $out/filtered.seqinfo.csv | sort) \
      <(grep '>' $out/filtered2.seqinfo.csv | sort)

    diff -q \
      <(grep '>' $out/filtered.details.csv | sort) \
      <(grep '>' $out/filtered2.details.csv | sort)

done
