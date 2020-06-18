FROM continuumio/miniconda3:4.8.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the luslab/group-nextflow-clip pipeline"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 procps=2:3.3.15-2 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install conda packages
COPY tools/bedtools/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-bedtools-0.1/bin:$PATH