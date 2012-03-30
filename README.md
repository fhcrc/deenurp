# De-novo reference set creation

Similarity-search based reference sequence selection

## Running

The `refset` package under the current directory provides to subcommands, accessed via the script `refset.py`:

### `refset.py search-sequences`

Searches a set of sequences against a one or more FASTA files containing possible reference sequences.

This subcommand does the following:

1. Clusters at a user specified identity
1. For each sequence database provided, searches all sequences (*not* cluster
   seeds) from clusters which don't have any hits. If any sequence within the
   cluster has a hit, sequences from the cluster are not searched against
   following databases.
1. Merges clusters which share a best_hit sequence
1. Saves the sequences and results in a sqlite database

### `refset.py select-references`

Given the output of `search-sequences`, `select-references` attempts to find an
optimal set of reference sequences.

1. For each sequence cluster, select a user-specified number of candidate
   references, preferring hits to highly weighted sequences.
1. Align the query sequences and reference candidates with cmalign
1. Use FastTree to make a tree of the reference sequences
1. Place the query sequences on the tree
1. Choose a specified number with `rppr voronoi`


# TODO
* Function to determine number of sequences to include for each cluster?
* Re-add type-strains (harder with some named, some unnamed sequences)
