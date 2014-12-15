#!/bin/bash

set -e

rm -rf output

DEENURP=${DEENURP-../../deenurp.py}
DATA=../../deenurp/test/data
mkdir -p output

for aligner in cmalign muscle vsearch; do
    time $DEENURP pairwise-distances --aligner $aligner \
	 $DATA/test_db_head.fasta output/${aligner}.csv
done
