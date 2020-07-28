#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

//Process definition
process homer_annotatePeaks {
    publishDir "${params.outdir}/homer",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-homer:latest'

    input:
        tuple val(sample_id), path(peaks)
        path(fasta)
        path(gtf)

    output:
        tuple val(sample_id), path("${sample_id}.annotatePeaks.txt")
    
    script:
    
    //SHELL
    """
    annotatePeaks.pl ${peaks} ${fasta} -gid -gtf ${gtf} > ${sample_id}.annotatePeaks.txt
    """
}

