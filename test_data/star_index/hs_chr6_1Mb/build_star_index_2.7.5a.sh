#!/bin/bash -ue

STAR \
  --runMode genomeGenerate \
  --runThreadN 2 \
  --genomeDir 2.7.5a \
  --genomeFastaFiles ../raw_genome/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa \
  --genomeSAindexNbases 9 \
  --outFileNamePrefix 2.7.5a