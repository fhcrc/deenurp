# `subcommands/expand_named.py`

## Purpose
Recruit additional named examples for underrepresented taxa using reciprocal
search against an unnamed sequence database.

## Inputs and Outputs
- Inputs:
  - named FASTA + seqinfo + taxonomy.
  - unnamed FASTA search database.
- Outputs:
  - `output.fasta` with appended accepted recruits.
  - `output.seq_info.csv` with inferred metadata rows.

## Main Flow
1. Find taxa below `--min-at-rank` at `--rank`.
2. Extract their current exemplar sequences.
3. Search exemplars against unnamed database.
4. Extract hits and search them back against exemplars.
5. Keep reciprocal-consistent hits and append to output datasets.

## Detailed Function Notes
- `uclust.search(...)`:
  - Runs `vsearch --usearch_global` with UC output.
  - Supports a broader `search_pct_id` then stricter post-filter cutoff.
- `uclust.parse_uclust_out(...)`:
  - Parses hit records used for forward and reverse matching logic.
- `util.ntf(...)`:
  - Manages temporary FASTA and UC intermediates.

## External Binaries
- `vsearch`.
