#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl = 2

// Process definition
process racon {
    label "max_cores"
    label "high_mem"
    label "regular_queue"
    
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/racon:1.4.12--he513fc3_0"

    input:
        val opts
        tuple val(meta), path(reads)
        path(overlap_paf)
        path(assembly_fasta)

    output:
        tuple val(meta), path("racon.fasta"), emit: fasta

    script:
	//Build the command line options
	racon_command = "racon --threads ${task.cpus} \
			$reads \
			$overlap_paf \
			$assembly_fasta > racon.fasta "

    //SHELL
    """
    ${racon_command}
    """
}
