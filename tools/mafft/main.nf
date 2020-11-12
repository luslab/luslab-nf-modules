#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process mafft {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/mafft:7.471--h516909a_0"

    input:
        val opts
        tuple val(meta), path (fasta)

    output:
        tuple val(meta), path ("*.mfa"), emit: mfa

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    mafft_command = "mafft $args --thread ${task.cpus} ${fasta} > ${fasta.simpleName}.mfa "

    if (params.verbose){
        println ("[MODULE] mafft command: " + mafft_command)
    }

    // SHELL
    """
    ${mafft_command}
    """
}
