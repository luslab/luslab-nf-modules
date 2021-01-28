#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process filtlong {
    label "low_cores"
    label "high_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/filtlong:0.1.0--0"

    input:
        val opts
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("*.filtlong.fastq"), emit: fastq

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    // Bear in mind that the output of filtlong is not compressed. You may want to do so using
    // another channel to not take up too much disk space.
    filtlong_command = "filtlong $args ${fastq} > ${fastq.simpleName}.filtlong.fastq"

    if (params.verbose){
        println ("[MODULE] filtlong command: " + filtlong_command)
    }
    //SHELL
    """
    ${filtlong_command}
    """
}
