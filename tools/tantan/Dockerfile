FROM nfcore/base:1.10
LABEL authors="charles.plessy@oist.jp" \
      description="Docker image containing all requirements for the nf-lra tantan module"

# Install apt packages
RUN apt-get update \
	&& apt-get install -y tantan \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Not installing from Bioconda, because tantan is version 13 in bioconda,
# but version 22 in Debian stable...
## Install conda packages
#COPY environment.yml /
#RUN conda env create -f /environment.yml && conda clean -a
#ENV PATH /opt/conda/envs/nfcore-module-flye/bin:$PATH
