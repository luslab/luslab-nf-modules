
FROM nfcore/base:1.7

LABEL authors="charlotte.west@crick.ac.uk" \
      description="Docker image containing requirements to run SEACR module"

# Update conda
RUN conda update -n base -c defaults conda

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-seacr/bin:$PATH