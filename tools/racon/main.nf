#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process racon {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-racon:latest"

    input:
        val opts
        tuple val(meta), path(reads), path(overlapPaf), path(assemblyFasta)

    output:
        tuple val(meta), path("racon.fasta"), emit: fasta

    script:
	//Build the command line options
	racon_command = "racon --threads ${task.cpus} \
			$reads \
			$overlapPaf \
			$assemblyFasta > racon.fasta " 

	//SHELL
    """
    ${racon_command}
    """
}
