#!/bin/bash

set -e

TESTS_DIR=$(dirname $BASH_SOURCE)

# to include a test, add a directory to the list below; lines
# beginning with '#' are skipped
while read subdir; do
    if echo $subdir | grep -qv -E '^#'; then
	echo $subdir
	(cd $TESTS_DIR/$subdir && ./run.sh)
	# (cd $TESTS_DIR/$subdir && bash -v ./run.sh)
    fi
done <<EOF
search-select
fill-lonely
pairwise-distances
hrefpkg-build
filter-outliers
EOF

exit
