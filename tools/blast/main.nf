#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process blast_makeblastdb {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*"), emit: blast_db

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    blast_command = "makeblastdb $args -in ${fasta}"

    if (params.verbose){
        println ("[MODULE] blast_makeblastdb command: " + minionqc_command)
    }

	//SHELL
    """
    ${blast_command}
    """
}
