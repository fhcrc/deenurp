FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt-get update && apt-get install --assume-yes \
    build-essential \
    gfortran \
    git \
    liblapack-dev \
    libopenblas-dev \
    make \
    python-dev \
    python-pip \
    python2.7 \
    unzip \
    wget

# Add files
RUN mkdir /usr/local/share/deenurp/
ADD bin /usr/local/share/deenurp/bin
ADD tests /usr/local/share/deenurp/tests
ADD deenurp /usr/local/share/deenurp/deenurp
ADD deenurp.py setup.py requirements.txt MANIFEST.in /usr/local/share/deenurp/

# Install deenurp and dependencies
RUN cd /usr/local/share/deenurp/ && \
  PYTHON=/usr/bin/python2.7 \
  DEENURP=/usr/local/share/deenurp/ \
  bin/bootstrap.sh /usr/local/

# Run tests
RUN cd /usr/local/share/deenurp && tests/run.sh
