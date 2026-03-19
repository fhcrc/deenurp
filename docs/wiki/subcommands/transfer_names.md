# `subcommands/transfer_names.py`

## Purpose
Transfer taxonomy names from an external FASTA+seqinfo set onto matching
reference-package sequences at or above a percent-identity threshold.

## Inputs and Outputs
- Inputs:
  - target reference package.
  - external FASTA + seqinfo + taxtable.
- Outputs:
  - updated `taxonomy` and `seq_info` resources in the refpkg.
  - optional audit log of rename decisions.

## Main Flow
1. Load new taxonomy and sequence metadata.
2. Load existing refpkg aligned sequences + metadata + taxonomy.
3. Search refpkg sequences against external FASTA.
4. For each high-identity hit, decide rename behavior on conflicts.
5. Add missing taxonomy nodes/ranks when transferring novel taxa.
6. Commit updated taxonomy/seq_info resources to refpkg.

## Detailed Function Notes
- `uclust.search(...)`:
  - Executes `vsearch --usearch_global` with requested identity.
- `uclust.parse_uclust_out(...)`:
  - Parses hit records used to drive rename operations.
- `add_to_taxonomy(...)`:
  - Walks up from missing node to nearest existing ancestor, then inserts
    missing nodes in parent->child order and updates rank list if needed.
- `util.as_fasta(...)` and `util.ntf(...)`:
  - Create temporary FASTA/CSV resources for search and refpkg updates.

## External Binaries
- `vsearch`.
