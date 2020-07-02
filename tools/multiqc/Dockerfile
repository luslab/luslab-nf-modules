FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for multiqc"

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/luslab-nf-modules-multiqc/bin:$PATH