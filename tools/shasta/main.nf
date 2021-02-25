#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process shasta {
    label "high_cores"
    label "max_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/shasta:0.6.0--hc9558a2_0"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("ShastaRun/Assembly.fasta"), emit: fasta
        tuple val(meta), path("ShastaRun/Assembly.gfa"), emit: gfa
        tuple val(meta), path("**{csv,dot,html,conf}"), emit: report

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    shasta_command = "shasta $args --threads ${task.cpus} --input ${reads}"

    //SHELL
    """
    ${shasta_command}
    """
}
