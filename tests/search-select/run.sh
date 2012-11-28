#!/bin/sh
set -u

RDP_BASE=../rdp_10_30_named1200bp_subset

deenurp search-sequences ground-reads.fa search.db $RDP_BASE.fasta $RDP_BASE.seqinfo.csv --group-field=tax_id
deenurp select-references search.db refs.fasta --seqinfo-out refs.seqinfo.csv --output-meta refs.meta.csv --mpi-args '-np 2'
