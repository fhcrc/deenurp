FROM ubuntu:20.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt-get update && apt-get install --assume-yes --no-install-recommends \
    build-essential \
    gfortran \
    git \
    liblapack-dev \
    libopenblas-dev \
    make \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    unzip \
    wget

# Add files
RUN mkdir /usr/local/share/deenurp/
ADD bin /usr/local/share/deenurp/bin
ADD tests /usr/local/share/deenurp/tests
ADD deenurp /usr/local/share/deenurp/deenurp
ADD deenurp.py setup.py requirements.txt /usr/local/share/deenurp/

# Install deenurp and dependencies
RUN cd /usr/local/share/deenurp/ && \
    PYTHON=/usr/bin/python3.8 \
    DEENURP=/usr/local/share/deenurp/ \
    bin/bootstrap.sh /usr/local/

# clean up sources apt packages
RUN rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /usr/local/share/deenurp/src && \
    apt-get purge -y --auto-remove \
      build-essential \
      unzip \
      git \
      python3-dev \
      make

# create some mount points
RUN mkdir -p /app /fh /mnt /run/shm

# Run tests
RUN python3 -m deenurp.test && \
    cd /usr/local/share/deenurp && \
    tests/run.sh
