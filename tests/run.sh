#!/bin/bash

for f in $(find . -wholename ./run.sh -prune -o -name run.sh -print); do
    echo $(dirname $f)
    (cd $(dirname $f) && ./run.sh)
done
