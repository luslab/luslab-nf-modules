#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process getcrosslinks {
    publishDir "${params.outdir}/get_crosslinks",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-get_crosslinks:latest'

    input:
      tuple val(sample_id), path(bam), path (fai)

    output:
      tuple val(sample_id), path ("${bam[0].simpleName}.xl.bed.gz"), emit: crosslinkBed

    script:

    //SHELL
    """
    bedtools bamtobed -i ${bam[0]} > dedupe.bed
    bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
    bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
    cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${bam[0].simpleName}.xl.bed.gz
    """
}