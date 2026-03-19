# `subcommands/cluster_refs.py`

## Purpose
Build clustered references and cluster labels by combining taxonomy-derived
named sequences with optional unnamed candidates.

## Inputs and Outputs
- Inputs:
  - `named_sequence_file` FASTA and matching `seqinfo`.
  - taxtable for taxonomy structure.
  - optional unnamed sequences and metadata.
- Outputs:
  - merged sequence FASTA (`sequence_out`).
  - sequence metadata with `cluster` column (`seqinfo_out`).

## Main Flow
1. Index input FASTA for fast random retrieval.
2. Load taxonomy and populate node sequence memberships from seqinfo.
3. Emit all sequences already at `--cluster-rank` as fixed clusters.
4. Build an intermediate pool of above-rank and unnamed sequences.
5. Remove redundant candidates close to already-named references.
6. Cluster remaining candidates into OTUs.
7. Write FASTA and seqinfo with per-sequence cluster assignment.

## Detailed Function Notes
- `wrap.read_seq_file(sequence_file)`:
  - Scans FASTA once and stores byte start/end offsets per sequence ID.
  - Used for random-access extraction without reparsing whole FASTA.
- `wrap.esl_sfetch(sequence_file, name_iter, output_fp, fa_idx)`:
  - Uses byte offsets from `fa_idx` to copy raw FASTA records.
  - This is pure Python I/O; it does not call the `esl-sfetch` binary.
- `uclust.search(...)`:
  - Calls `vsearch --usearch_global` and writes UC-format search results.
  - Supports `maxaccepts`, `maxrejects`, and post-filter modes.
- `uclust.cluster(...)`:
  - Calls `vsearch --cluster_fast` or `--cluster_smallmem`.
  - Produces UC cluster output for seed/hit relationships.
- `uclust.parse_uclust_out(...)` and `uclust.sequences_by_cluster(...)`:
  - Parse UC lines and group records by cluster seed.

## External Binaries
- `vsearch` (search + clustering).
