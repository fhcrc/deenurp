#!/bin/bash
set -u
set -e

RDP_BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

rm -rf output
mkdir -p output

$DEENURP fill_lonely --exclude-taxids exclude_taxids.txt --include-taxids include_taxids.txt $RDP_BASE.fasta $RDP_BASE.seqinfo.csv $RDP_BASE.taxonomy.csv refs.fasta refs.seqinfo.csv output/fill_lonely.fa output/fill_lonely.csv
