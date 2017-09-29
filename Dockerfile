FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt-get update && apt-get install -y \
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
RUN mkdir /bin/deenurp
ADD . /bin/deenurp

# Install deenurp and dependencies
RUN cd /bin/deenurp/bin/ && \
    chmod +x bootstrap.sh && \
    PYTHON=/usr/bin/python2.7 DEENURP=/bin/deenurp /bin/bash bootstrap.sh
RUN cd /bin/deenurp/bin && /bin/bash -c "source bin-env/bin/activate"

# Add to path
ENV PATH="/bin/deenurp/bin/bin-env/bin:${PATH}"

# Run tests
RUN cd /bin/deenurp/tests && /bin/bash run.sh
