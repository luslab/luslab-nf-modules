FROM nfcore/base:1.7
LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core velocyto module"

# Update and install system libraries
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends  \
    # zlib1g-dev \
    # libbz2-dev \
    # liblzma-dev \
    gcc \
    # mono-mcs \
    build-essential \
    # cmake \
    # less \
    # libbamtools-dev \
    # libboost-dev \
    # libboost-iostreams-dev \
    # libboost-log-dev \
    # libboost-system-dev \
    # libboost-test-dev \
    # libcurl4-openssl-dev \
    # libssl-dev \
    # libz-dev \
    # nano \
    pandoc && \
    rm -rf /var/lib/apt/lists/*
    
# Update conda
RUN conda update -n base -c defaults conda

# Set channels
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge

# Install conda packages
RUN conda install -y python=3.6 \
    numpy=1.19.1 \
    scipy=1.2.0 \
    cython=0.29.4 \
    numba=0.41.0 \
    matplotlib=3.0.2 \
    scikit-learn=0.20.2 \
    h5py=2.9.0 \
    click=7.0

RUN pip install pysam==0.15.2 \
    velocyto==0.17.17
