#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process htseq_count {
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    tag "${meta.sample_id}"
    
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", overwrite: true
    
    container "quay.io/biocontainers/htseq:0.12.4--py37h70f9b12_1"

    input:
        val opts
        tuple val(meta), path(reads)
        path gtf

    output:
        tuple val(meta), path("${prefix}"), emit: counts

    script:
        // Check main args string exists and strip whitespace
        args = ''
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        // Construct CL line
        htseq_command = "htseq-count ${opts.args} $reads $gtf > ${prefix}"

        // Log
        if (params.verbose){
            println ("[MODULE] htseq_count command: " + htseq_command)
        }

        //SHELL
        """
        ${htseq_command}
        """
}