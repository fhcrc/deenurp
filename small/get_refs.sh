#!/bin/sh

set -e
set -u

../refset.py search-sequences dedup.fasta small.db --weights dedup_info.csv ../../../../micro_refset/corpus/build/rdp.fasta
../refset.py select-references small.db test_refs.fasta
