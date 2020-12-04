#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process purge_haplotigs {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "quay.io/biocontainers/purge_haplotigs:1.1.1--0"
    container "luslab/nf-modules-purge_haplotigs:1.0.0"

    input:
        val opts
        tuple val(meta), path(bam)
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*"), emit: whatever

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    purge_haplotigs_command = "purge_haplotigs $args readhist -b ${bam} -g ${fasta}"

    //SHELL
    """
    ${purge_haplotigs_command}
    """
}
