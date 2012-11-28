#!/bin/sh
BASE=rdp_10_30_named1200bp_subset
mkdir -p hrefpkg
deenurp hrefpkg-build --index-rank=family $BASE.fasta $BASE.seqinfo.csv $BASE.taxonomy.csv --output-dir hrefpkg
