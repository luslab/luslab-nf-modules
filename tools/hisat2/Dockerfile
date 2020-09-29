FROM nfcore/base:1.7
LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core hisat2 module"

COPY VERSION .

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-hisat2/bin:$PATH