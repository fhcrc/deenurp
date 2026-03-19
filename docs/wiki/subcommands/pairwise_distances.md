# `subcommands/pairwise_distances.py`

## Purpose
Generate a labeled pairwise distance matrix for input sequences.

## Inputs and Outputs
- Input: FASTA file (`seqs`).
- Output: CSV distance matrix (`distmat`) with sequence IDs in header.

## Main Flow
1. Derive output prefix from FASTA filename.
2. Dispatch by `--aligner`:
   - `vsearch`: direct pairwise identity pipeline.
   - `cmalign`: CM alignment then FastTree matrix.
   - `muscle`: MSA then FastTree matrix.
3. Save matrix with `numpy.savetxt`.

## Detailed Function Notes
- `filter_outliers.distmat_pairwise(...)`:
  - Uses `wrap.vsearch_allpairs_files(...)` and blast6 parsing.
- `filter_outliers.distmat_cmalign(...)`:
  - Uses `wrap.cmalign_files(...)`, Stockholm->FASTA conversion,
    and `outliers.fasttree_dists(...)`.
- `filter_outliers.distmat_muscle(...)`:
  - Uses `wrap.muscle_files(...)` then `outliers.fasttree_dists(...)`.

## External Binaries
- `vsearch`
- `cmalign` + `FastTree`
- `muscle` + `FastTree`
