#!/bin/bash

# create a singularity image

if [[ -z $1 ]]; then
    echo "usage: $0 <TAG>"
    docker images deenurp
    exit
fi

tag="$1"
img="deenurp:$tag"

docker run \
       -v "${DOCKER_HOST#unix://}:/var/run/docker.sock" \
       -v "$(pwd):/output" \
       --privileged -t --rm \
       quay.io/singularity/docker2singularity \
       --name "deenurp-${tag}.sif" "$img"


