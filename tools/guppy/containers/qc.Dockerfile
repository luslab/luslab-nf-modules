FROM nfcore/base:1.7
LABEL authors="candice.hermant@crick.ac.uk" \
      description="Docker image containing all requirements for guppy quality control process"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-guppy-qc/bin:$PATH