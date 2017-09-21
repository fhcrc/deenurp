FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt-get update && apt-get install -y build-essential wget unzip make gfortran libopenblas-dev liblapack-dev python2.7 python-dev git python-pip

# Add files
RUN mkdir /bin/deenurp
ADD . /bin/deenurp

# Install conda and hdbscan
ENV PATH="/root/miniconda2/bin/:${PATH}"
RUN cd /bin/deenurp && \
	wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
	bash Miniconda2-latest-Linux-x86_64.sh -b && \
	conda install -c conda-forge hdbscan


# Install python requirements
RUN pip install -r /bin/deenurp/requirements.txt

# Install deenurp
RUN cd /bin/deenurp/bin/ && chmod +x bootstrap.sh && PYTHON=/usr/bin/python2.7 DEENURP=/bin/deenurp ./bootstrap.sh
RUN cd /bin/deenurp/bin/ && /bin/bash -c "source bin-env/bin/activate" 

# Add to path
ENV PATH="/bin/deenurp/bin/bin-env/bin/:${PATH}"

# Run tests
RUN cd /bin/deenurp/tests && bash run.sh