#!/bin/bash
set -u
set -e

RDP_BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

rm -rf output
mkdir -p output

$DEENURP search_sequences ground-reads.fa output/search.db \
	 $RDP_BASE.fasta \
	 $RDP_BASE.seqinfo.csv \
	 --group-field=tax_id \
	 --blacklist=blacklist.txt
$DEENURP select_references output/search.db output/refs.fasta \
	 --seqinfo-out output/refs.seqinfo.csv \
	 --output-meta output/refs.meta.csv \
	 --min-mass-prop 0.01 \
	 --whitelist whitelist.txt
