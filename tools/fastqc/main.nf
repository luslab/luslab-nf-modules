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
    mv_command = ''

    // Check the report name
    if(params.fastqc_reportname && params.fastqc_reportname != ''){
        mv_command = 
        "mv ${reads.simpleName}*.html ${sample_id}_${params.fastqc_reportname}.html && mv ${reads.simpleName}*.zip ${sample_id}_${params.fastqc_reportname}.zip"
    }

    // Log
    if (params.verbose){
        println ("[MODULE] fastqc command: " + fastqc_command)
    }
    
    //SHELL
    """
    ${fastqc_command} 
    ${mv_command}
    """
}