#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

//Process definition
process bedtools_intersect {
    publishDir "${params.outdir}/bedtools/intersect",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-bedtools:latest'

    input: 
        tuple val(sample_id), path(reads), path(regions_file)

    output: 
        tuple val(sample_id), path("${sample_id}.annotated.bed"), emit: annotatedBed

    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.bedtools_intersect_args && params.bedtools_intersect_args != '') {
        ext_args = params.bedtools_intersect_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    intersect_command = "bedtools intersect -a ${regions_file} -b $reads ${args} > ${sample_id}.annotated.bed"

    // Log
    if (params.verbose){
        println ("[MODULE] bedtools/intersect command: " + intersect_command)
    }

    //SHELL
    """
    ${intersect_command}
    """
}