#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process flye {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-flye:latest"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_id}/assembly.fasta"), emit: assemblyFasta

    script:
	//Build the command line options
    flye_command = "flye --genome-size ${opts.genome_size} \
			--threads ${task.cpus} \
			--out-dir ${meta.sample_id} \
			--nano-raw $reads"

	//SHELL
    """
    ${flye_command}
    """
}
