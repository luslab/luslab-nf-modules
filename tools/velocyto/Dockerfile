FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core bedtools module"

RUN apt-get update && apt-get install -y \
    software-properties-common

RUN add-apt-repository ppa:deadsnakes/ppa

RUN apt install \
    python3.7 \
    python3-pip


# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-velocyto/bin:$PATH

RUN pip install pysam
RUN pip install velocyto

# # Install conda packages
# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a
# ENV PATH /opt/conda/envs/nfcore-module-cutadapt/bin:$PATH