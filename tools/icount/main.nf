#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process icount {
    publishDir "${params.outdir}/icount",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-icount:latest'

    input:
        tuple val(sample_id), path(bed), path(seg)

    output:
        tuple val(sample_id), path("${bed.simpleName}.xl.peaks.bed.gz"), emit: peaks
        tuple val(sample_id), path("${bed.simpleName}.scores.tsv"), emit: peak_scores
        tuple val(sample_id), path("${bed.simpleName}.xl.clusters.bed.gz"), emit: clusters
    
    script:

    //SHELL
    """
    iCount peaks $seg $bed ${bed.simpleName}.xl.peaks.bed.gz \
        --scores ${bed.simpleName}.scores.tsv \
        --half_window ${params.icount_half_window} \
        --fdr ${params.icount_fdr}

    zcat ${bed.simpleName}.xl.peaks.bed.gz | \
    bedtools merge -i stdin -s -d ${params.icount_half_window} -c 4,5,6 -o distinct,sum,distinct | \
    gzip > ${bed.simpleName}.xl.clusters.bed.gz
    """
}