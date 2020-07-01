FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk,sviatoslav.sidorov@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core STAR module"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-umitools/bin:$PATH
