# `subcommands/deduplicate_sequences.py`

## Purpose
Collapse identical sequence strings into representatives while preserving and
aggregating metadata.

## Inputs and Outputs
- Inputs:
  - sequence FASTA.
  - seqinfo CSV indexed by `seqname`.
- Outputs:
  - deduplicated FASTA.
  - representative seqinfo with computed `weight`.

## Main Flow
1. Compute SHA1 hash for each sequence string.
2. Join `seqhash` into seqinfo.
3. Group by `seqhash` and optional `--group-by` keys.
4. Choose representative row per group (`tail(1)` after optional sorting).
5. Set `weight = group_size` and write outputs.

## Detailed Function Notes
- `util.file_opener(mode=...)`:
  - Opens plain, `.gz`, or `.bz2` files with a common interface.
- `util.Counter(iterable, ...)`:
  - Wraps iterables and prints progress counters to stderr.

## External Binaries
None.
