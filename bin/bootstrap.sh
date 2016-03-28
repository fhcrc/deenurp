#!/bin/bash

set -e

# Usage: [PYTHON=path/to/python] [DEENURP=path/to/deenurp] bootstrap.sh
#
# Create a virtualenv, and install requirements to it.
#
# specify a python interpreter using
# `PYTHON=path/to/python bootstrap.sh`
# specify path to the deenurp source directory using
# `DEENURP=path/to/deenurp bootstrap.sh`

# installs deenurp and dependencies to $VIRTUAL_ENV if defined

# Will attempt to install python packages from wheels if $PIP_FIND_LINKS is defined
# and pip --use-wheel is specified

# Will attempt to create wheels if $PIP_WHEEL_DIR is defined
# see https://pip.pypa.io/en/latest/user_guide.html#environment-variables

# use values of $PIP_WHEEL_DIR and $PIP_FIND_LINKS from environment if
# possible, otherwise use a location in src
SRCDIR=$(readlink -f src)
export PIP_WHEEL_DIR=${PIP_WHEEL_DIR-$SRCDIR/cache/pip/wheels}
export PIP_FIND_LINKS=${PIP_FIND_LINKS-file://$PIP_WHEEL_DIR}

mkdir -p src
mkdir -p $PIP_WHEEL_DIR

srcdir(){
    tar -tf $1 | head -1
}

if [[ -n "$1" ]]; then
    venv="$1"
elif [[ -n $VIRTUAL_ENV ]]; then
    venv=$VIRTUAL_ENV
else
    venv=$(basename $(pwd))-env
fi

if [[ -z $PYTHON ]]; then
    PYTHON=$(which python)
fi

# Defines the default source directory for deenurp as the parent of
# the directory containing this file.
if [[ -z $DEENURP ]]; then
    DEENURP=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
fi

VENV_VERSION=1.11.6
PPLACER_BUILD=1.1.alpha17
INFERNAL_VERSION=1.1.1
UCLUST_VERSION=1.2.22
RAXML_VERSION=8.0.5
MUSCLE_VERSION=3.8.31
VSEARCH_VERSION=1.9.5

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

# full path; set by activate
venv=$VIRTUAL_ENV

# install pplacer and accompanying python scripts
PPLACER_DIR=pplacer-Linux-v${PPLACER_BUILD}
PPLACER_ZIP=${PPLACER_DIR}.zip

pplacer_is_installed(){
    $venv/bin/pplacer --version 2> /dev/null | grep -q "$PPLACER_BUILD"
}

if pplacer_is_installed; then
    echo -n "pplacer is already installed: "
    $venv/bin/pplacer --version
else
    mkdir -p src && \
	(cd src && \
	wget -nc https://github.com/matsen/pplacer/releases/download/v$PPLACER_BUILD/$PPLACER_ZIP && \
	unzip -o $PPLACER_ZIP && \
	cp $PPLACER_DIR/{pplacer,guppy,rppr} $venv/bin && \
	pip install -U $PPLACER_DIR/scripts)
    # confirm that we have installed the requested build
    if ! pplacer_is_installed; then
	echo -n "Error: you requested pplacer build $PPLACER_BUILD "
	echo "but $($venv/bin/pplacer --version) was installed."
	echo "Try removing src/$PPLACER_ZIP first."
	exit 1
    fi
fi

# install infernal and easel binaries
INFERNAL=infernal-${INFERNAL_VERSION}-linux-intel-gcc

if [ ! -f $venv/bin/cmalign ]; then
    (cd src && \
	wget -nc http://eddylab.org/infernal/${INFERNAL}.tar.gz && \
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
	wget -nc http://drive5.com/uclust/uclustq${UCLUST_VERSION}_i86linux64 && \
	chmod +x uclustq${UCLUST_VERSION}_i86linux64 && \
	ln -f uclustq${UCLUST_VERSION}_i86linux64 uclust)
else
    echo "$(uclust --version) is already installed"
fi

# install FastTree and FastTreeMP
if [ ! -f $venv/bin/FastTree ] | [ ! -f $venv/bin/FastTreeMP ]; then
    (cd $venv/bin && \
	wget -nc http://www.microbesonline.org/fasttree/FastTree && \
	wget -nc http://www.microbesonline.org/fasttree/FastTreeMP && \
	chmod +x FastTree*)
else
    echo "FastTree is already installed: $(FastTree -expert 2>&1 | head -1)"
fi

# install raxmlHPC-SSE3 and raxmlHPC-PTHREADS-SSE3
if [ ! -f $venv/bin/raxmlHPC-SSE3 ] | [ ! -f $venv/bin/raxmlHPC-PTHREADS-SSE3 ]; then
    (cd src && \
	wget -nc https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz && \
	tar -xf v${RAXML_VERSION}.tar.gz && \
	cd standard-RAxML-${RAXML_VERSION} && \
	rm -f *.o && make -f Makefile.SSE3.gcc && \
	rm -f *.o && make -f Makefile.SSE3.PTHREADS.gcc && \
	mv raxmlHPC-SSE3 raxmlHPC-PTHREADS-SSE3 $venv/bin
    )
else
    echo "raxml is already installed: $(raxmlHPC-SSE3 | grep RAxML)"
fi

echo $venv

# install VSEARCH
vsearch_is_installed(){
    $venv/bin/vsearch --version | grep -q "$VSEARCH_VERSION"
}

if vsearch_is_installed; then
    echo -n "vsearch is already installed: "
    $venv/bin/vsearch --version
else
    (cd src && \
	    wget -nc https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
	    tar -xf vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
	    cp vsearch-${VSEARCH_VERSION}-linux-x86_64/bin/vsearch $venv/bin && \
	    chmod +x $venv/bin/vsearch
    )
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
	    wget -nc http://www.drive5.com/muscle/downloads${MUSCLE_VERSION}/muscle${MUSCLE_VERSION}_src.tar.gz && \
	    tar -xf muscle${MUSCLE_VERSION}_src.tar.gz && \
	    cd muscle${MUSCLE_VERSION}/src && \
	    ./mk && cp muscle $venv/bin)
fi

# install wheels library
# if [ -n "$PIP_WHEEL_DIR" ]; then
#   pip install wheel
#   pip wheel wheel
# fi

# required to build and cache wheels - doing this greatly speeds up travis CI tests
pip install -U pip
pip install -U wheel
pip wheel wheel

# install python requirements; note that `pip install -r
# requirements.txt` fails due to install-time dependencies.
while read line; do
  pip wheel "$line" # --verbose helps with the travis test timeouts
  pip install --upgrade "$line"
done < "$DEENURP/requirements.txt"

pip install -e "$DEENURP"

# correct any more shebang lines
virtualenv --relocatable $venv
