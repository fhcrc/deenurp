#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(clstutils))
suppressPackageStartupMessages(library(ape))

# Find outliers from an alignment
# Arguments: <alignment_file> <cutoff> <prune_out>
args <- commandArgs(TRUE)
alignment <- args[1]
cutoff <- as.numeric(args[2])
prune_out <- args[3]

a <- read.dna(alignment, format='fasta')
dm <- dist.dna(a, pairwise.deletion=TRUE, as.matrix=TRUE)
prune <- findOutliers(dm, cutoff=cutoff)
to.prune <- colnames(dm)[prune]
cat(paste(to.prune, collapse='\n'), file=prune_out)
