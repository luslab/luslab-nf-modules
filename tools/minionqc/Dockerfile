FROM nfcore/base:1.7

LABEL authors="michael.mansfield@oist.jp" \
      description="Docker image containing all requirements for the nf-core minionqc module"

COPY VERSION .

# Install all of R's dependencies, then install R
RUN apt-get update \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-minionqc/bin:$PATH
