# `subcommands/hrefpkg_build.py`

## Purpose
Construct a hierarchy of per-taxon reference packages and an index package
that points to child packages.

## Inputs and Outputs
- Inputs: full sequence FASTA, seqinfo CSV, taxonomy taxtable.
- Outputs:
  - multiple `<tax_id>.refpkg` directories.
  - `index.refpkg` and `index.csv` mapping tax_id->refpkg path.
  - train/test FASTA and not-in-package FASTA side products.

## Main Flow
1. Load taxonomy and seqinfo.
2. Optionally partition taxonomy for validation mode.
3. For each index-rank node, build a child refpkg in parallel.
4. Build `index.refpkg` by merging child seqinfo/taxonomy scope.
5. Extract sequences not represented in any generated package.

## Detailed Function Notes
- `tax_id_refpkg(...)` internals:
  - choose sequences per sub-taxon (`PER_TAXON`),
  - align via `wrap.cmalign(...)`,
  - tree via `wrap.fasttree(...)`,
  - write refpkg resources and reroot package.
- `wrap.esl_sfetch(...)`:
  - Extracts selected IDs into FASTA intermediates.
- `util.require_executable('rppr')`:
  - Enforces availability before `Refpkg.reroot()` workflow.

## External Binaries
- `cmalign`
- `FastTree` / `FastTreeMP`
- `rppr`
