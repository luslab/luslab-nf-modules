#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process last {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

    input:
        val opts
        tuple val(meta), path(reads),
        path fasta_file

    output:
        tuple val(meta), path("*.maf"), emit: mappedReads

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = ""
        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}