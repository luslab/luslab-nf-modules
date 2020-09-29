#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process nanoplot {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-nanoplot:latest"

    input:
        val opts
        tuple val(meta), path(nanopore_reads)

    output:
        tuple val(meta), path("*{pdf,html,log}"), emit: nanoplotOutputs

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    nanoplot_command = "NanoPlot -t ${task.cpus} ${args} -p nanoplot. --fastq $nanopore_reads"

    if (params.verbose){
        println ("[MODULE] nanoplot command: " + nanoplot_command)
    }

	//SHELL
    """
    ${nanoplot_command}
    """
}
