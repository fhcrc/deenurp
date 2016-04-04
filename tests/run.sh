#!/bin/bash

set -e

vsearch --version

(cd tests/pairwise-distances && ./run.sh)
(cd tests/search-select && ./run.sh)
exit

TESTS_DIR=$(dirname $BASH_SOURCE)

for f in $(find $TESTS_DIR -wholename $TESTS_DIR/run.sh -prune -o -name run.sh -print); do
    echo $(dirname $f)
    (cd $(dirname $f) && ./run.sh)
done
