#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Main process
process multiqc {
    publishDir "${params.outdir}/multiqc",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-multiqc:latest'
    //container 'quay.io/biocontainers/multiqc:1.9--py_1'

    input:
      path(reports)

    output:
      path "multiqc_report.html", emit: report
      path "multiqc_data/multiqc.log", emit: log
        
    script:

    // Check for custom args
    args = ""
    if(params.multiqc_args && params.multiqc_args != '') {
        ext_args = params.multiqc_args
        args += " " + ext_args.trim()
    }

    // Check for custom config
    custom_config_file = ''
    if(params.multiqc_config && params.multiqc_config != '') {
        $custom_config_file = "--config ${params.multiqc_config}"
    }

    multiqc_command = "multiqc . $args $custom_config_file"

    // Log
    if (params.verbose){
        println ("[MODULE] multiqc command: " + multiqc_command)
    }

    """
    ${multiqc_command}
    """
}