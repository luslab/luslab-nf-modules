#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Trimming reusable component
process cutadapt {
    publishDir "${params.outdir}/cutadapt",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-cutadapt:latest'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fq.gz"), emit: trimmedReads
        path "*.cutadapt.txt", emit: report

    script:

    // Init
    args = "--log=${sample_id}.cutadapt.txt"

    // Check main args string exists and strip whitespace
    if(params.cutadapt_args) {
        ext_args = params.cutadapt_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    cutadapt_command = "cutadapt ${args} $reads -o ${sample_id}.trimmed.fq.gz > ${sample_id}_cutadapt.txt"

    // Log
    if (params.verbose){
        println ("[MODULE] cutadapt command: " + cutadapt_command)
    }

    //SHELL
    """
    ${cutadapt_command} 
    """
}
