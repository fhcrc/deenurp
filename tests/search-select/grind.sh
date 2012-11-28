#!/bin/sh

grinder -rf ../hrefpkg-build/rdp_10_30_named1200bp_subset.fasta -fr primer.fasta \
  -length_bias 0 \
  -unidirectional 1 \
  -copy_bias 0 \
  -tr 1000 \
  -bn ground
