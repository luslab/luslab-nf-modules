#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process def
process bowtie2_align {
    publishDir "${params.outdir}/bowtie2",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        tuple val(sample_id), path(reads), path(index)

    //output:
    //    tuple val(sample_id), path("*.trimmed.fq.gz"), emit: trimmedReads
    //    path "*.txt", emit: report

    script:
        // Check main args string exists and strip whitespace
        args = ''
        if(params.bowtie2_args && params.bowtie2_args != '') {
            ext_args = params.bowtie2_args
            args += " " + ext_args.trim()
        }

        """
        bowtie2
        """
}