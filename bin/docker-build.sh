#!/bin/bash

set -e

if ! git diff-index --quiet HEAD; then
    echo "Error: git repo is dirty - commit and try again"
    echo
    git status
    exit 1
fi

# remove any lingering pyc files
find . -name '*.pyc' -exec rm {} \;

# docker tag should reflect current version number
tag=deenurp:$(git describe --tags --dirty)
docker build -t $tag .
