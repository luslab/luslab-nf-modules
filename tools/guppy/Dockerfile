FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04
LABEL author="Candice Hermant" \
      description="Docker image containing all requirements for the guppy module"

# 
RUN --runtime=nvidia --name luslab/nf-modules-guppy -i -t -v <local_dir_of_fast5:/media/>  \
    nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04 /bin/bash

# Get the guppy dev file
RUN wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_3.1.5-1~xenial_amd64.deb

# Install guppy packages
RUN apt-get install libbost-all-dev \
    apt-get remove ont-guppy

# Install guppy dev file
RUN dpkg -i --ignore-depends=nvidia-384,libcuda1-384 ont_guppy_3.1.5-1~xenial_amd64.deb
