#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process augustus {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/augustus:3.3.3--pl526h0faeac2_5"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*"), emit: augustus

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    augustus_command = "augustus $args --outfile=${fasta.simpleName}.augustus --species=${opts.species} ${fasta}"

    if (params.verbose){
        println ("[MODULE] augustus command: " + augustus_command)
    }

    //SHELL
    """
    ${augustus_command}
    """
}
