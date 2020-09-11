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
				//you'll need to change these
        //tuple val(meta), path(reads)

    output:
				// you'll also need to change these
        //tuple val(meta), path("${meta.sample_id}/assembly.fasta"), emit: assemblyFasta
				//tuple path("${meta.sample_id}/assembly_info.txt"), path("${meta.sample_id}/flye.log"), emit: log

    script:
	//Build the command line options
	//change these
  //  racon_command = "flye --genome-size ${opts.genome_size} \
	//		--threads ${task.cpus} \
	//		--out-dir ${meta.sample_id} \
	//		--nano-raw $reads"

	//SHELL
    """
    ${racon_command}
    """
}
