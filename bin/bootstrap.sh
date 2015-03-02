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
VSEARCH_VERSION=1.1.1

check_version(){
    # usage: check_version module version-string
    "$PYTHON" <<EOF 2> /dev/null
import $1
from distutils.version import LooseVersion
assert LooseVersion($1.__version__) >= LooseVersion("$2")
EOF
}

# create virtualenv if necessary, downloading source if available
# version is not up to date.
VENV_URL="https://pypi.python.org/packages/source/v/virtualenv"
if [[ ! -f "${venv:?}/bin/activate" ]]; then
    # if the system virtualenv is up to date, use it
    if check_version virtualenv $VENV_VERSION; then
	echo "using $(which virtualenv) (version $(virtualenv --version))"
    	virtualenv "$venv"
    else
	echo "downloading virtualenv version $VENV_VERSION"
	if [[ ! -f src/virtualenv-${VENV_VERSION}/virtualenv.py ]]; then
	    mkdir -p src
	    (cd src && \
		wget -N ${VENV_URL}/virtualenv-${VENV_VERSION}.tar.gz && \
		tar -xf virtualenv-${VENV_VERSION}.tar.gz)
	fi
	"$PYTHON" src/virtualenv-${VENV_VERSION}/virtualenv.py "$venv"
    fi
else
    echo "virtualenv $venv already exists"
fi

source $venv/bin/activate
# make relocatable (necessary to prevent long shebang lines)
virtualenv --relocatable $venv

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

# install VSEARCH
vsearch_is_installed(){
    $venv/bin/vsearch --version | grep -q "$VSEARCH_VERSION"
}

if vsearch_is_installed; then
    echo -n "vsearch is already installed: "
    $venv/bin/vsearch --version
else
    (cd src && \
	    wget -N https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-linux-x86_64 && \
	    mv vsearch-${VSEARCH_VERSION}-linux-x86_64 $venv/bin && \
	    chmod +x $venv/bin/vsearch-${VSEARCH_VERSION}-linux-x86_64 && \
	    ln -f $venv/bin/vsearch-${VSEARCH_VERSION}-linux-x86_64 $venv/bin/vsearch)
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
virtualenv --relocatable $venv
