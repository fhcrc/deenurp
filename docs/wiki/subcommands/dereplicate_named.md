# `subcommands/dereplicate_named.py`

## Purpose
Dereplicate named references by clustering sequences within each grouping key
(default `species`) and keeping cluster seeds.

## Inputs and Outputs
- Inputs: named FASTA, seqinfo CSV, optional taxonomy taxtable.
- Outputs:
  - representative FASTA (`--seqs-out`).
  - optional derep map (`--derep-map-out`).
  - optional filtered seqinfo (`--seq-info-out`).

## Main Flow
1. Load seqinfo and optional taxonomy columns.
2. Group rows by `--group-on`.
3. For singleton groups, mark the single sequence as seed.
4. For larger groups, extract FASTA subset and run clustering at `--id`.
5. Combine seeds across groups and write requested outputs.

## Detailed Function Notes
- `wrap.read_seq_file(...)`:
  - Creates FASTA offset index for efficient repeated fetches.
- `wrap.esl_sfetch(...)`:
  - Extracts representative sequences from the source FASTA.
- `uclust.cluster(...)`:
  - Runs `vsearch` clustering and writes UC records.
- `uclust.parse_uclust_as_df(...)`:
  - Parses UC into DataFrame; normalizes seed target labels.

## External Binaries
- `vsearch`.
