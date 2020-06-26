FROM continuumio/miniconda3:4.8.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core fastqc module"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-fastqc/bin:$PATH