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

    // Check main args string exists and strip whitespace
    args = ""
    if(params.fastqc_args && params.fastqc_args != '') {
        ext_args = params.fastqc_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    fastqc_command = "fastqc${args} --threads ${task.cpus} $reads"


    // Log
    if (params.verbose){
        println ("[MODULE] fastqc command: " + fastqc_command)
    }
    
    //SHELL
    readList = reads.collect{it.toString()}
    if(readList.size > 1){
        if(params.fastqc_reportname && params.fastqc_reportname != ''){
            """
            ${fastqc_command}
            mv ${reads[0].simpleName}*.html ${reads[0].simpleName}_${params.fastqc_reportname}.html
            mv ${reads[0].simpleName}*.zip ${reads[0].simpleName}_${params.fastqc_reportname}.zip
            mv ${reads[1].simpleName}*.html ${reads[1].simpleName}_${params.fastqc_reportname}.html
            mv ${reads[1].simpleName}*.zip ${reads[1].simpleName}_${params.fastqc_reportname}.zip
            """
        } else {
            """
            ${fastqc_command}
            """
        }
    } else {
        if(params.fastqc_reportname && params.fastqc_reportname != ''){
            """
            ${fastqc_command}
            mv ${reads[0].simpleName}*.html ${reads[0].simpleName}_${params.fastqc_reportname}.html
            mv ${reads[0].simpleName}*.zip ${reads[0].simpleName}_${params.fastqc_reportname}.zip
            """
        } else {
            """
            ${fastqc_command}
            """
        }
    }
}