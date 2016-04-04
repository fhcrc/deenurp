#!/bin/bash

set -e

TESTS_DIR=$(dirname $BASH_SOURCE)

while read subdir; do
    if echo $subdir | grep -qv -E '^#'; then
	echo $subdir
	(cd $TESTS_DIR/$subdir && ./run.sh)
    fi
done <<EOF
search-select
fill-lonely
pairwise-distances
hrefpkg-build
filter-outliers
EOF

exit

TESTS_DIR=$(dirname $BASH_SOURCE)

for f in $(find $TESTS_DIR -wholename $TESTS_DIR/run.sh -prune -o -name run.sh -print); do
    echo $(dirname $f)
    (cd $(dirname $f) && ./run.sh)
done
