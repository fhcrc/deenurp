FROM python:3.11-slim-bookworm
LABEL org.opencontainers.image.authors="sminot@fredhutch.org,nhoffman@uw.edu,crosenth@uw.edu"

# Install prerequisites
RUN apt-get update && \
apt-get upgrade --assume-yes && \
apt-get install --assume-yes --no-install-recommends build-essential git unzip wget

# Add files
RUN mkdir /usr/local/share/deenurp/
COPY bin /usr/local/share/deenurp/bin
COPY tests /usr/local/share/deenurp/tests
COPY deenurp /usr/local/share/deenurp/deenurp
COPY deenurp.py setup.py requirements.txt /usr/local/share/deenurp/

# Install deenurp and dependencies
RUN cd /usr/local/share/deenurp/ && \
    PYTHON=/usr/local/bin/python3 \
    DEENURP=/usr/local/share/deenurp/ \
    bin/bootstrap.sh /usr/local/ \
    pip install --upgrade --requirement requirements.txt

# clean up sources apt packages
RUN rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /usr/local/share/deenurp/src && \
    apt-get purge -y --auto-remove git

# create some mount points
RUN mkdir -p /app /fh /mnt /run/shm
