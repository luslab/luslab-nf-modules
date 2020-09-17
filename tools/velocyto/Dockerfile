FROM quay.io/biocontainers/samtools:1.10--h2e538c0_3

FROM quay.io/biocontainers/velocyto.py:0.17.17--py37h97743b1_2

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements for the nf-core velocyto module"
      
# copy samtools from first biocontainer
COPY --from=0 /usr/local/bin/samtools /bin/