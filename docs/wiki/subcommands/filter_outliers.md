# `subcommands/filter_outliers.py`

## Purpose
Filter likely outlier/misannotated references within each taxonomic group
using pairwise distances and configurable outlier rules.

## Inputs and Outputs
- Inputs:
  - sequence FASTA, seqinfo CSV, taxonomy taxtable.
  - optional previous detailed results.
- Outputs:
  - required filtered FASTA (`--output-seqs`).
  - optional filtered seqinfo and detailed seqinfo with diagnostics.

## Main Flow
1. Load taxonomy and assign sequence IDs to taxonomic nodes.
2. Keep sequences classified above filter rank.
3. For each taxon at filter rank:
   - skip/filter trivially for rare taxa or `--no-filter` cases, or
   - extract taxon FASTA and compute distance matrix.
4. Mark outliers via radius or cluster strategy.
5. Merge outcomes across taxa and emit filtered outputs.

## Detailed Function Notes
- `wrap.vsearch_allpairs_files(...)`:
  - Calls `vsearch --allpairs_global ... --blast6out ...`.
  - Produces blast6-like pairwise identity table.
- `wrap.cmalign_files(...)` + `outliers.fasttree_dists(...)`:
  - `cmalign` generates alignment; `FastTree -makematrix` converts to dists.
- `wrap.muscle_files(...)` + `outliers.fasttree_dists(...)`:
  - Alternative alignment + same FastTree matrix generation.
- `outliers.outliers(...)`:
  - Medoid-distance thresholding.
- `outliers.outliers_by_cluster(...)`:
  - Cluster-aware pruning using SciPy/HDBSCAN clustering state.
- `wrap.esl_sfetch(...)`:
  - Fetches per-taxon FASTA subsets and final kept set.

## External Binaries
- `vsearch` (pairwise mode)
- `cmalign` (cmalign mode)
- `muscle` (muscle mode)
- `FastTree` (matrix from alignment)
