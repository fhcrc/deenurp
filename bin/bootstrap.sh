#!/bin/bash

# Usage: [PYTHON=path/to/python] [DEENURP=path/to/deenurp] bootstrap.sh [virtualenv-name]
#
# Create a virtualenv, and install requirements to it.
#
# specify a python interpreter using
# `PYTHON=path/to/python bootstrap.sh`
# specify path to the deenurp source directory using
# `DEENURP=path/to/deenurp bootstrap.sh`

set -e

srcdir(){
    tar -tf $1 | head -1
}

if [[ -z $1 ]]; then
    venv=$(basename $(pwd))-env
else
    venv=$1
fi

if [[ -z $PYTHON ]]; then
    PYTHON=$(which python)
fi

if [[ -z $DEENURP ]]; then
    DEENURP="."
fi

mkdir -p src

VENV_VERSION=1.11.2
PPLACER_VERSION=1.1
INFERNAL_VERSION=1.1
UCLUST_VERSION=1.2.22

# Create a virtualenv using a specified version of the virtualenv
# source. This also provides setuptools and pip. Inspired by
# http://eli.thegreenplace.net/2013/04/20/bootstrapping-virtualenv/
VENV_URL='http://pypi.python.org/packages/source/v/virtualenv'

# download virtualenv source if necessary
if [ ! -f src/virtualenv-${VENV_VERSION}/virtualenv.py ]; then
    (cd src && \
	wget -N ${VENV_URL}/virtualenv-${VENV_VERSION}.tar.gz && \
	tar -xf virtualenv-${VENV_VERSION}.tar.gz)
fi

# create virtualenv if necessary
if [ ! -f $venv/bin/activate ]; then
    $PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py $venv
    $PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py --relocatable $venv
else
    echo "found existing virtualenv $venv"
fi

source $venv/bin/activate

# full path; set by activate
venv=$VIRTUAL_ENV

# install python requirements; note that `pip install -r
# requirements.txt` fails due to install-time dependencies.
while read line; do
    pip install -U "$line"
done < "$DEENURP/requirements.txt"

# install deenurp
pip install -e "$DEENURP"

# install pplacer and accompanying python scripts
PPLACER_TGZ=pplacer-v${PPLACER_VERSION}-Linux.tar.gz
if [ ! -f $venv/bin/pplacer ]; then
	(cd src && \
	wget -N http://matsen.fhcrc.org/pplacer/builds/$PPLACER_TGZ && \
	tar -xf $PPLACER_TGZ && \
	cp $(srcdir $PPLACER_TGZ)/{pplacer,guppy,rppr} $venv/bin && \
	pip install -U $(srcdir $PPLACER_TGZ)/scripts && \
	rm -r $(srcdir $PPLACER_TGZ))
else
      echo "$(pplacer --version) is already installed"
fi

# install infernal and easel binaries
INFERNAL=infernal-${INFERNAL_VERSION}-linux-intel-gcc
venv_abspath=$(readlink -f $venv)

if [ ! -f $venv/bin/cmalign ]; then
    (cd src && \
	wget -N http://selab.janelia.org/software/infernal/${INFERNAL}.tar.gz && \
	for binary in cmalign cmconvert esl-alimerge esl-sfetch; do
	    tar xvf ${INFERNAL}.tar.gz --no-anchored binaries/$binary
	done && \
	    cp ${INFERNAL}/binaries/* $venv/bin && \
	    rm -r ${INFERNAL}
    )
else
    echo "cmalign is already installed: $(cmalign -h | sed -n 2p)"
fi

# install uclust
if [ ! -f $venv/bin/uclust ]; then
    (cd $venv/bin && \
	wget -N http://drive5.com/uclust/uclustq${UCLUST_VERSION}_i86linux64 && \
	chmod +x uclustq${UCLUST_VERSION}_i86linux64 && \
	ln -f uclustq${UCLUST_VERSION}_i86linux64 uclust)
else
    echo "$(uclust --version) is already installed"
fi

# install FastTree and FastTreeMP
if [ ! -f $venv/bin/FastTree ] | [ ! -f $venv/bin/FastTreeMP ]; then
    (cd $venv/bin && \
	wget -N http://www.microbesonline.org/fasttree/FastTree && \
	wget -N http://www.microbesonline.org/fasttree/FastTreeMP && \
	chmod +x FastTree*)
else
    echo "FastTree is already installed: $(FastTree -expert 2>&1 | head -1)"
fi

# correct any more shebang lines
$PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py --relocatable $venv
