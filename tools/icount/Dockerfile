FROM continuumio/miniconda3:4.8.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core icount module"
# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2 \
 unzip=6.0-23+deb10u1 \
 procps=2:3.3.15-2 \
 build-essential=12.6 \
 zlib1g-dev=1:1.2.11.dfsg-1 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-icount-0.1/bin:$PATH

# Install iCount
WORKDIR /home
RUN wget https://github.com/tomazc/iCount/archive/master.zip && mv master.zip icount.zip && unzip icount.zip
WORKDIR /home/iCount-master
# hadolint ignore=SC2102,DL3013
RUN pip install -e .[test]