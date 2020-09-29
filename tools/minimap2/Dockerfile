FROM nfcore/base:1.7
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core minimap2 module"

COPY VERSION .

# Change dir
WORKDIR /home

# Install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN mv /home/minimap2-2.17_x64-linux /home/minimap2

# Add to path
ENV PATH /home/minimap2:$PATH
