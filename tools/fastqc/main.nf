#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process fastqc {
    publishDir "${params.outdir}/fastqc",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-fastqc:latest'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path ("*.{zip,html}"), emit: fastqcOutput
        path "*.{zip,html}", emit: report

    script:

    // Init
    args = "--log=${sample_id}.{zip,html}"

    // Check main args string exists and strip whitespace
    if(params.fastqc_args) {
        ext_args = params.fastqc_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    fastqc_command = "fastqc ${args} --threads ${task.cpus} $reads"

    // Log
    if (params.verbose){
        println ("[MODULE] fastqc command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    """
}