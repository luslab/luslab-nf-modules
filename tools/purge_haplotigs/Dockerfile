FROM nfcore/base:1.7

LABEL authors="michael.mansfield@oist.jp" \
      description="Docker image containing all requirements for the nf-core purge_haplotigs module"

COPY VERSION .

RUN apt-get update \
    && apt-get install -y --no-install-recommends libpng-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-purge_haplotigs/bin:$PATH

# Add python script to the container
COPY purge_haplotigs_minima.py /

WORKDIR /work
