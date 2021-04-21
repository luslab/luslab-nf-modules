#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process porechop {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/porechop:0.2.4--py38hed8969a_1"

    input:
        val opts
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("*.porechop.fastq.gz"), emit: fastq

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    porechop_command = "porechop $args --threads ${task.cpus} -i ${fastq} -o ${fastq.simpleName}.porechop.fastq.gz"

    //SHELL
    """
    ${porechop_command}
    """
}
