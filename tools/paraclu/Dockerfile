FROM continuumio/miniconda3:4.8.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core paraclu module"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2 \
 unzip=6.0-23+deb10u1 \
 procps=2:3.3.15-2 \
 build-essential=12.6 \
 zlib1g-dev=1:1.2.11.dfsg-1 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-paraclu-0.1/bin:$PATH

# Install Paraclu
WORKDIR /home
RUN mkdir bin && wget http://cbrc3.cbrc.jp/~martin/paraclu/paraclu-9.zip && unzip paraclu-9.zip
WORKDIR /home/paraclu-9
RUN make && cp paraclu /home/bin/paraclu && cp paraclu-cut.sh /home/bin/paraclu-cut.sh
