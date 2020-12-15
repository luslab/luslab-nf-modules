FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core peka module"

COPY VERSION .

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 apt-utils=1.4.11  \
 build-essential=12.3 \
 zlib1g-dev=1:1.2.8.dfsg-5 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Update conda
RUN conda update -n base -c defaults conda

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-peka/bin:$PATH

# Install PEKA
WORKDIR /home
RUN git clone https://github.com/ulelab/imaps.git --branch 0.5.1 --single-branch
WORKDIR /home/imaps
RUN pip install -e .
WORKDIR /home
RUN mkdir src && mv imaps/src/imaps/sandbox/kmers.py src/kmers.py
