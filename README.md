# Installation

First, install binary dependencies:

* uclust
* easel tools
* cmalign
* pplacer
* tax2tree (and pycogent)
* BioPython
* FastTree
* dnaclust (switch to uclust?)

# De-novo reference set creation

Similarity-search based reference sequence selection

## Running

The `refset` package under the current directory provides to subcommands,
accessed via the script `refset.py`, or the command `refset` if installed.

Subcommands fall into two general categories:

* Building a set of reference sequences for use in refpkg building
* Selecting sequences for a specific reference package

## Creating a sequence set for refpkg building

### `refset filter-outliers`

Removes sequences from a reference database that are more than a specified
distance from the centroid of their tax id.

### `refset expand-named`

Expands poorly-represented names in a sequence file by similarity search

### `refset cluster-refs`

Cluster reference sequences, first by tax-id at a specified rank (default:
species), then by similarity for unnamed sequences or sequences not classified
to the desired rank.  Serves as input to `search-sequences`.

## Selecting sequences for a reference package

### `refset hrefpkg-build`

Builds a set of hierarchical reference packages.

### `refset search-sequences`

Searches a set of sequences against a FASTA file containing possible reference sequences.

This subcommand does searches sequences against a reference FASTA file, saving
the results and some metadata to a sqlite database for use in
`select-references`

### `refset select-references`

Given the output of `search-sequences`, `select-references` attempts to find a
good set of reference sequences.

For each reference cluster  with a minimal amount of sequences having best hits
to the cluster, (see `cluster-refs`), selects a set number of sequences to
serve as references.

### `refset add-reps`

Fetches sequences from a sequence file which match the taxtable for a reference
set at a given rank. Useful for adding type strains.

### `refset tax2tree`

Runs the `tax2tree` program on a reference package, updating the `seq_info`
file.

Sequences whose lineage changes are relabeled. The prior `tax_id` is added to
the `seq_info` file in the reference package.

# TODO
* Fix loneliness

