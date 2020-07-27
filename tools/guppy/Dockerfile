# Set the base image to Ubuntu 16.04
FROM ubuntu:16.04

# File Author / Maintainer
LABEL author="Charlotte Berthelier <bertheli@biologie.ens.fr>" 

ARG PACKAGE_VERSION=3.6.0
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes $BUILD_PACKAGES && \
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_cpu_${PACKAGE_VERSION}-1~xenial_amd64.deb && \
    apt-get install -y --no-install-recommends libzmq5=4.1.4-7ubuntu0.1 libhdf5-cpp-11=1.8.16+docs-4ubuntu1.1 libcurl4-openssl-dev=7.47.0-1ubuntu2.15 libssl-dev=1.0.2g-1ubuntu4.16 libhdf5-10=1.8.16+docs-4ubuntu1.1 libboost-regex1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-log1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-atomic1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-chrono1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-date-time1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-filesystem1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-program-options1.58.0=1.58.0+dfsg-5ubuntu3.1 libboost-iostreams1.58.0=1.58.0+dfsg-5ubuntu3.1 && \
    dpkg -i *.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*