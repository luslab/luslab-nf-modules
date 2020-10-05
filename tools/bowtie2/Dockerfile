#FROM nfcore/base:1.7

FROM quay.io/biocontainers/bioconvert:0.4.3--py_0

FROM quay.io/biocontainers/samtools:1.10--h2e538c0_3

FROM quay.io/biocontainers/bowtie2:2.4.1--py38h1c8e9b9_3

LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf bowtie2 module"

COPY --from=0 /usr/local/bin/pigz /usr/local/bin

COPY --from=1 /usr/local/bin/samtools /usr/local/bin

COPY VERSION .

# Install conda packages
#COPY environment.yml /
#RUN conda env create -f /environment.yml && conda clean -a
#ENV PATH /opt/conda/envs/nfcore-module-bowtie2/bin:$PATH