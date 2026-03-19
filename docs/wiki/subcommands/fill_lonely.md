# `subcommands/fill_lonely.py`

## Purpose
Add additional representatives for taxa that are alone under their parent at
a chosen rank.

## Inputs and Outputs
- Inputs:
  - full search FASTA + seqinfo + taxtable.
  - chosen FASTA + seqinfo to augment.
- Outputs:
  - augmented FASTA.
  - augmented seqinfo rows for added representatives.

## Main Flow
1. Identify lonely nodes (`--lonely-rank`) in chosen taxonomy.
2. For each lonely node, gather sibling-taxon sequences from full taxonomy.
3. Align sibling sequences, build tree, and prune to `--number-of-reps`.
4. Union all additional reps (plus include/exclude taxid overrides).
5. Append FASTA and seqinfo rows.

## Detailed Function Notes
- `wrap.read_seq_file(...)` and `wrap.esl_sfetch(...)`:
  - Build/read offset index and fetch selected records.
- `wrap.cmalign(sequences, ...)`:
  - Runs `cmalign` alignment and yields aligned SeqRecords.
- `wrap.fasttree(...)`:
  - Runs `FastTreeMP` when available for multithreading, otherwise `FastTree`.
- `wrap.rppr_min_adcl_tree(...)`:
  - Executes `rppr min_adcl_tree` and returns leaf names to prune.
- `util.ntf(...)`:
  - Temporary storage for FASTA/tree intermediates.

## External Binaries
- `cmalign`
- `FastTree` or `FastTreeMP`
- `rppr`
