Bootstrap:docker
From:ubuntu:18.04

%setup
  mkdir ${SINGULARITY_ROOTFS}/src

%files
  bin /src/bin
  tests /src/tests
  deenurp /src/deenurp
  setup.py /src/
  requirements.txt /src/
  MANIFEST.in /src/

%post
  apt-get update && apt-get install --assume-yes --no-install-recommends \
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
  cd /src/ && \
    PYTHON=/usr/bin/python2.7 \
    DEENURP=/src/ \
    bin/bootstrap.sh /usr/local/
  rm -rf /var/lib/apt/lists/*
  rm -rf /root/.cache/pip
  rm -rf /src/src
  apt-get purge -y --auto-remove \
    build-essential \
    unzip \
    git \
    python-dev \
    make
  mkdir -p /fh /app /mnt  # create bind points
