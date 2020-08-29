FROM nfcore/base:1.7
LABEL authors="anob.chakrabarti@crick.ac.uk" \
      description="Docker image containing all requirements for the luslab-nf-modules STAR module"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/luslab-nf-star/bin:$PATH