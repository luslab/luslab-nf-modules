FROM continuumio/miniconda
LABEL authors="candice.hermant@crick.ac.uk" \
      description="Docker image containing all requirements for the Minimap2 module"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-minimap2/bin:$PATH
