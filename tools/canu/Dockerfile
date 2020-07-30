FROM nfcore/base:1.7
LABEL authors="Luscombe Lab" \
      description="Docker image containing all requirements for the nf-core canu module"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 build-essential=12.3 \
 zlib1g-dev=1:1.2.8.dfsg-5 \
 && apt-get clean -y \
 && rm -rf /var/lib/apt/lists/*

# Download canu
RUN wget https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Linux-amd64.tar.xz
RUN tar -xJf canu-2.0.*.tar.xz
RUN ln -s canu-2.0/*/bin/canu /canu
