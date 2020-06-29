#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Trimming reusable component
process samtools_index {
    publishDir "${params.outdir}/samtools/index",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*.bam.bai"), emit: baiFiles
 
    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.samtools_index_args && params.samtools_index_args != '') {
        ext_args = params.samtools_index_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    index_command = "samtools index ${args} -@ ${task.cpus} $bam"

    // Log
    if (params.verbose){
        println ("[MODULE] samtools/index command: " + index_command)
    }
    
    """
    ${index_command}
    """
}
