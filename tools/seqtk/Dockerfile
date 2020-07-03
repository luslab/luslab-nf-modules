FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core seqtk module"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 build-essential=12.3 \
 zlib1g-dev=1:1.2.8.dfsg-5 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Setup exec folder
WORKDIR /home
RUN mkdir bin
ENV PATH /home/bin:$PATH

# Install seqtk
RUN git clone https://github.com/lh3/seqtk.git
WORKDIR /home/seqtk
RUN make && cp seqtk /home/bin/seqtk && chmod 755 seqtk

