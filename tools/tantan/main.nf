#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl = 2

// Process definition
process tantan {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-tantan:latest"

    input:
        val opts
        tuple val(meta), path(genome)

    output:
        tuple val(meta), path("${meta.sample_id}.${opts.suffix}"), emit: tantanRepeats

    script:
	//Build the command line options
    command = "tantan -w${opts.max_period} -f4 ${genome} > ${meta.sample_id}.${opts.suffix}"
	//SHELL
    """
    ${command}
    """
}
