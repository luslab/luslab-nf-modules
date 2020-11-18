#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process emboss_seqret {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/emboss:6.6.0--h847b9ba_3"

    input:
        val opts
        tuple val(meta), path(in_seq)

    output:
        tuple val(meta), path ("*"), emit: out_seq

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    emboss_command = "seqret $args -sequence ${in_seq} -sformat1 ${opts.input_format} -outseq ${in_seq.simpleName}.${opts.output_format} -osformat2 ${opts.output_format}"

    if (params.verbose){
        println ("[MODULE] EMBOSS command: " + emboss_command)
    }

    // SHELL
    """
    ${emboss_command}
    """
}
