#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process multiqc {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/multiqc:1.9--py_1'

    input:
        val(opts)
        path(multiqc_config)
        path(reports)

    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data
        
    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        multiqc_command = "multiqc -f $args ."
        if(opts.custom_config) {
            multiqc_command = "multiqc -f$args -c $multiqc_config ."
        }

        if (params.verbose){
            println ("[MODULE] multiqc command: " + multiqc_command)
        }

        """
        ${multiqc_command}
        """
}