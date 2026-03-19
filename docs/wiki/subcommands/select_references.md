# `subcommands/select_references.py`

## Purpose
Choose final reference sequence set from the search DB, balancing cluster
coverage and query-supported abundance.

## Inputs and Outputs
- Inputs:
  - search DB from `search-sequences`.
  - include/exclude cluster or sequence lists (optional).
- Outputs:
  - selected FASTA.
  - optional selection metadata and subset seqinfo.

## Main Flow
1. Parse inclusion/exclusion options.
2. Work from a temporary copy of DB to avoid in-place mutation issues.
3. Load indexed FASTA offsets for query/reference files.
4. Iterate selected refs from `select.choose_references(...)`.
5. Deduplicate by ID and by exact sequence string.
6. Write FASTA and optional metadata artifacts.

## Detailed Function Notes
- `select.choose_references(...)`:
  - Reads hit weights per cluster/sample from DB.
  - Fetches cluster reference and query members.
  - Submits per-cluster tasks to select up to `--refs-per-cluster` refs.
- `select.select_sequences_for_cluster(...)`:
  - Pre-clusters candidate refs using `uclust.cluster` (`vsearch`).
  - Aligns refs+queries with `wrap.cmalign`.
  - Builds temporary refpkg with `wrap.as_refpkg` (FastTree-backed).
  - Places reads with `wrap.pplacer`, applies `wrap.guppy_redup`, and prunes
    with `wrap.rppr_min_adcl`.
- `wrap.esl_sfetch(...)`:
  - Fetches cluster member FASTA records by ID.

## External Binaries
- `vsearch`
- `cmalign`
- `FastTree` / `FastTreeMP`
- `pplacer`
- `guppy`
- `rppr`
