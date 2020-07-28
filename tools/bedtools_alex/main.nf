#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

//Process definition
process bedtools_intersect {
    publishDir "${params.outdir}/bedtools/intersect",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-bedtools:latest'

    input:
        tuple val(sample_id1), path(reads1)
        tuple val(sample_id2), path(reads2)

    output:
        path("intersected.bed")

    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.bedtools_intersect_args && params.bedtools_intersect_args != '') {
        ext_args = params.bedtools_intersect_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    intersect_command = "bedtools intersect -a ${reads1} -b ${reads2} ${args} > intersected.bed"

    // Log
    if (params.verbose){
        println ("[MODULE] bedtools/intersect command: " + intersect_command)
    }

    //SHELL
    """
    ${intersect_command}
    """
}