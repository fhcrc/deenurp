# `subcommands/cmalign.py`

## Purpose
Run Infernal covariance-model alignment and return aligned FASTA plus optional
score table.

## Inputs and Outputs
- Input: unaligned FASTA file.
- Outputs:
  - aligned FASTA file.
  - optional CSV of cmalign scores (`--scores`).

## Main Flow
1. Allocate temporary Stockholm output.
2. Call `wrap.cmalign_files(...)` to run `cmalign`.
3. Convert Stockholm alignment to FASTA.
4. Return score DataFrame; optionally write it to CSV.

## Detailed Function Notes
- `wrap.cmalign_files(input_file, output_file, cm=CM, cpu=...)`:
  - Verifies `cmalign` is installed and version is Infernal 1.1.
  - Executes `cmalign --noprob --dnaout --cpu ... -o ...`.
  - Parses stdout score text into a typed pandas DataFrame.
- `util.ntf(...)`:
  - Temporary-file context manager that unlinks file on exit.

## External Binaries
- `cmalign`.
