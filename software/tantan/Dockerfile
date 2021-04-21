FROM nfcore/base:1.10
LABEL authors="charles.plessy@oist.jp" \
      description="Docker image containing all requirements for the nf-lra tantan module"

COPY VERSION .

# Install apt packages
RUN apt-get update \
	&& apt-get install --no-install-recommends -y tantan=22-1 \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Not installing from Bioconda, because tantan is version 13 in bioconda,
# but version 22 in Debian stable, on which nfcore/base:1.10 is built.
## Install conda packages
#COPY environment.yml /
#RUN conda env create -f /environment.yml && conda clean -a
#ENV PATH /opt/conda/envs/nfcore-module-flye/bin:$PATH
