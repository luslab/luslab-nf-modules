FROM nfcore/base:1.7
LABEL authors="michael.mansfield@oist.jp" \
      description="Docker image containing all requirements for the nf-lra flye module"

COPY VERSION .

# Install apt packages
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-flye/bin:$PATH
