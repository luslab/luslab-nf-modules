#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process hmmer_hmmscan {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/hmmer:3.3--he1b5a44_1"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.hmmscan"), emit: hmmscan
        tuple val(meta), path("*.tbl"), emit: table

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    hmmscan_command = "hmmscan $args --cpu ${task.cpus} --domtblout ${fasta.simpleName}.tbl ${opts.db} ${fasta} > ${fasta.simpleName}.hmmscan"

    if (params.verbose){
        println ("[MODULE] hmmscan_command command: " + hmmscan_command)
    }
    //SHELL
    """
    ${hmmscan_command}
    """
}

process hmmer_hmmsearch {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/hmmer:3.3--he1b5a44_1"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.hmmsearch"), emit: hmmsearch
        tuple val(meta), path("*.tbl"), emit: table

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    hmmer_hmmsearch = "hmmscan $args --cpu ${task.cpus} --domtblout ${fasta.simpleName}.tbl ${opts.db} ${fasta} > ${fasta.simpleName}.hmmsearch"

    if (params.verbose){
        println ("[MODULE] hmmer_hmmsearch command: " + hmmer_hmmsearch)
    }
    //SHELL
    """
    ${hmmer_hmmsearch}
    """
}
