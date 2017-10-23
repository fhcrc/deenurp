#!/bin/bash

# build from a local docker image

echo 'this is nor working for singularity 2.4'
exit 1

set -e

if [[ -z $1 ]]; then
    echo "usage: ./build.sh TAG [<outdir>]"
    echo
    echo "choose from one of the following tags:"
    echo
    docker image list deenurp
    exit 1
fi

outdir=$(readlink -f ${2-.})

version=$1
docker_image=deenurp:$1
docker_name=$(echo $docker_image | tr : _)_$RANDOM
img=$(readlink -f $outdir/deenurp-${version}.img)

if [[ -f $img ]]; then
    echo "$img already exists"
    exit 1
fi

# supposedly how it should be done for singularity 2.4
# must disable user namespaces with --userns=host to use --privileged
docker run \
-v /var/run/docker.sock:/var/run/docker.sock \
-v $outdir:/output \
--userns=host --privileged -t --rm \
singularityware/docker2singularity \
$docker_image

# singularity create --size 2300 $img
# docker run --name $docker_name $docker_image /bin/true
# docker export $docker_name | sudo singularity import $img
# docker rm $docker_name
# sudo singularity exec --writable $img mkdir -p /fh /app /mnt
