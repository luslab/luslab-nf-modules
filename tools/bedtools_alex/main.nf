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
        tuple val("${sample_id1}+${sample_id2}"), path("${sample_id1}+${sample_id2}_intersected.bed")
    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.bedtools_intersect_args && params.bedtools_intersect_args != '') {
        ext_args = params.bedtools_intersect_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    intersect_command = "bedtools intersect -a ${reads1} -b ${reads2} ${args} > '${sample_id1}+${sample_id2}_intersected.bed'"

    // Log
    if (params.verbose){
        println ("[MODULE] bedtools/intersect command: " + intersect_command)
    }

    //SHELL
    """
    ${intersect_command}
    """
}


//Process definition
process bedtools_subtract {
    publishDir "${params.outdir}/bedtools/subtract",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-bedtools:latest'

    input:
        tuple val(sample_id1), path(reads1)
        tuple val(sample_id2), path(reads2)

    output:
        tuple val("${sample_id1}-${sample_id2}"), path("${sample_id1}-${sample_id2}_subtracted.bed")

    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.bedtools_subtract_args && params.bedtools_subtract_args != '') {
        ext_args = params.bedtools_subtract_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    subtract_command = "bedtools subtract -a ${reads1} -b ${reads2} ${args} > '${sample_id1}-${sample_id2}_subtracted.bed'"

    // Log
    if (params.verbose){
        println ("[MODULE] bedtools/subtract command: " + subtract_command)
    }

    //SHELL
    """
    ${subtract_command}
    """
}