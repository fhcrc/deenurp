#!/bin/bash

# Usage: [PYTHON=path/to/python] [DEENURP=path/to/deenurp] bootstrap.sh [virtualenv-name]
#
# Create a virtualenv, and install requirements to it.
#
# specify a python interpreter using
# `PYTHON=path/to/python bootstrap.sh`
# specify path to the deenurp source directory using
# `DEENURP=path/to/deenurp bootstrap.sh`

# Will attempt to install python packages from wheels (using `pip
# install --use-wheel --find-links="$WHEELHOUSE"`) if a directory or
# url is defined in the environment variable `$WHEELHOUSE`

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

# Defines the default source directory for deenurp as the parent of
# the directory containing this file.
if [[ -z $DEENURP ]]; then
    DEENURP=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
fi

mkdir -p src

VENV_VERSION=1.11.6
PPLACER_BINARY_VERSION=1.1
PPLACER_BUILD=1.1.alpha16
INFERNAL_VERSION=1.1
UCLUST_VERSION=1.2.22
RAXML_VERSION=8.0.5
MUSCLE_VERSION=3.8.31

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

# install pplacer and accompanying python scripts
PPLACER_TGZ=pplacer-v${PPLACER_BINARY_VERSION}-Linux.tar.gz

pplacer_is_installed(){
    $venv/bin/pplacer --version 2> /dev/null | grep -q "$PPLACER_BUILD"
}

if pplacer_is_installed; then
    echo -n "pplacer is already installed: "
    $venv/bin/pplacer --version
else
    mkdir -p src && \
	(cd src && \
	wget -N http://matsen.fhcrc.org/pplacer/builds/$PPLACER_TGZ && \
	tar -xf $PPLACER_TGZ && \
	cp $(srcdir $PPLACER_TGZ)/{pplacer,guppy,rppr} $venv/bin && \
	pip install -U $(srcdir $PPLACER_TGZ)/scripts)
    # confirm that we have installed the requested build
    if ! pplacer_is_installed; then
	echo -n "Error: you requested pplacer build $PPLACER_BUILD "
	echo "but $($venv/bin/pplacer --version) was installed."
	echo "Try removing src/$PPLACER_TGZ first."
	exit 1
    fi
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

# install raxmlHPC-SSE3 and raxmlHPC-PTHREADS-SSE3
if [ ! -f $venv/bin/raxmlHPC-SSE3 ] | [ ! -f $venv/bin/raxmlHPC-PTHREADS-SSE3 ]; then
    (cd src && \
	wget -N https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz && \
	tar -xf v${RAXML_VERSION}.tar.gz && \
	cd standard-RAxML-${RAXML_VERSION} && \
	rm -f *.o && make -f Makefile.SSE3.gcc && \
	rm -f *.o && make -f Makefile.SSE3.PTHREADS.gcc && \
	mv raxmlHPC-SSE3 raxmlHPC-PTHREADS-SSE3 $venv/bin
    )
else
    echo "raxml is already installed: $(raxmlHPC-SSE3 | grep RAxML)"
fi

# install MUSCLE
muscle_is_installed(){
    $venv/bin/muscle -version | grep -q "$MUSCLE_VERSION"
}

if muscle_is_installed; then
    echo -n "muscle is already installed: "
    $venv/bin/muscle -version
else
    (cd src && \
	    wget -N http://www.drive5.com/muscle/downloads${MUSCLE_VERSION}/muscle${MUSCLE_VERSION}_src.tar.gz && \
	    tar -xf muscle${MUSCLE_VERSION}_src.tar.gz && \
	    cd muscle${MUSCLE_VERSION}/src && \
	    ./mk && cp muscle $venv/bin)
fi

# install python requirements; note that `pip install -r
# requirements.txt` fails due to install-time dependencies.
while read line; do
    if [[ -z "$WHEELHOUSE" ]]; then
	pip install -U "$line"
    else
	pip install --use-wheel --find-links="$WHEELHOUSE" "$line"
    fi
done < "$DEENURP/requirements.txt"

# install deenurp
pip install -e "$DEENURP"

# correct any more shebang lines
$PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py --relocatable $venv
