#FROM nfcore/base:1.7
FROM debian:unstable
LABEL authors="michael.mansfield@oist.jp" \
      description="Docker image containing all requirements for the nf-core minionqc module"

# Install all of R's dependencies, then install R
RUN apt-get update \
  && apt-get install -y libblas-dev libcurl4-openssl-dev libicu-dev liblapack-dev libpcre2-dev libxml2-dev libssl-dev wget curl python3-pip ruby ruby-dev \
  && apt-get install -y r-base \
	&& apt-get install -y r-cran-ggplot2

# Install R packages
RUN R -e "install.packages('data.table', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('futile.logger', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('plyr', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('readr', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('reshape2', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('scales', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('viridis', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('yaml', dependencies=TRUE, repos='https://cloud.r-project.org/')"

# Get the actual minionQC code
RUN curl https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R > /MinIONQC.R && chmod +x /MinIONQC.R
