# `subcommands/rdp_sequence_filter.py`

## Purpose
Deprecated command for splitting sequences into named/unnamed groups after
simple quality filters.

## Inputs and Outputs
- Inputs: FASTA file and seqinfo CSV in matching order.
- Outputs: named and unnamed FASTA + seqinfo files.

## Main Flow
1. Iterate FASTA records and seqinfo rows in lockstep.
2. Compute ambiguous-base proportion and length.
3. Reject records failing thresholds.
4. Route accepted records to named or unnamed outputs by
   `taxid_classified` flag.

## Detailed Function Notes
- `util.file_opener(...)`: opens compressed/uncompressed streams.
- `util.Counter(...)`: progress indicator wrapper for long iterations.

## External Binaries
None.
