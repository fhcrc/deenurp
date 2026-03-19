# `subcommands/add_reps.py`

## Purpose
Augment a sequence corpus by selecting FASTA records whose taxonomy intersects
with tax IDs observed in a taxtable at a target rank.

## Inputs and Outputs
- Inputs:
  - FASTA file to pull from.
  - `tax_map` CSV mapping `seqid,tax_id`.
  - taxonomy SQLite DB (`taxit new_database` output).
  - taxtable CSV.
- Output:
  - FASTA written to `outfile` containing selected representatives.

## Main Flow
1. Parse taxtable rows and collect tax IDs where `rank == --rank`.
2. Load `tax_map` and validate unique assignment per sequence.
3. For each mapped sequence tax ID, resolve lineage through taxonomy DB.
4. Keep sequence IDs whose lineage value at `--rank` is in selected tax IDs.
5. Stream FASTA and emit only selected sequences.

## Detailed Function Notes
- `wrap.load_tax_maps(fps, has_header=False)`:
  - Reads one or more tax-map CSV files.
  - Returns a dict of `{seqname: tax_id}`.
  - Raises `ValueError` if the same `seqname` maps to conflicting tax IDs.

## External Binaries
None in `deenurp` helper calls used by this module.
