#!/bin/bash
set -u
set -e

RDP_BASE=../rdp_10_30_named1200bp_subset
DEENURP=${DEENURP-../../deenurp.py}

$DEENURP search-sequences ground-reads.fa search.db $RDP_BASE.fasta $RDP_BASE.seqinfo.csv --group-field=tax_id
$DEENURP select-references search.db refs.fasta --seqinfo-out refs.seqinfo.csv --output-meta refs.meta.csv --mpi-args '-np 2'
