# `subcommands/search_sequences.py`

## Purpose
Build a SQLite search database linking query sequences to candidate reference
hits, including hit-ranking metadata used by later selection steps.

## Inputs and Outputs
- Inputs:
  - query FASTA (`sequence_file`).
  - reference FASTA (`ref_database`).
  - reference metadata CSV (`ref_meta`) with grouping field.
  - optional weight/sample map files and blacklist.
- Output:
  - SQLite DB (`output`) with schema + sequences + best-hit records.

## Main Flow
1. Parse optional `--sample-map`, `--weights`, and `--blacklist`.
2. Open SQLite output connection.
3. Call `search.create_database(...)` with thresholds and options.

## Detailed Function Notes
- `search.load_sample_map(fp)`:
  - Reads two-column query->sample mapping.
- `search.dedup_info_to_counts(fp, sample_map)`:
  - Converts guppy dedup rows into nested weight dict
    `{seqid: {sample: weight}}`.
- `search.create_database(...)` does the heavy lifting:
  1. Validates threshold relationship (`search_threshold` cannot exceed
     `search_identity` beyond tiny float tolerance).
  2. Calls `_create_tables(...)`:
     - executes SQL schema script from `deenurp/data/search.schema`,
     - writes run params (`fasta_file`, `ref_fasta`, `ref_meta`, identity,
       maxaccepts/maxrejects, group field, etc.) into `params` table.
  3. Calls `_load_sequences(...)`:
     - parses query FASTA,
     - inserts each sequence into `sequences`,
     - inserts sample-weight rows into `sequences_samples`.
  4. Calls `_search(...)`:
     - loads reference group membership from `ref_meta`,
     - runs similarity search via `uclust.search(...)`,
     - filters UC hits by configured identity,
     - applies `select_hits(...)` window around top hit,
     - enforces one kept hit per cluster,
     - inserts kept rows into `best_hits` and `ref_seqs`.

## External Binaries
- `vsearch` (indirectly via `uclust.search(...)` inside `_search(...)`).
