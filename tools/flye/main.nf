#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process flye {
    publishDir "${params.outdir}/flye",
        mode: "copy", overwrite: true

    container "luslab/nf-modules-flye:latest"

    input:
        tuple val(sample_id), val(genome_size), path(reads)

    output:
        tuple val(sample_id), path("_assembly"), emit: flyeAssembly

    script:

		//Build the command line options
    flye_command = "flye --genome-size ${genome_size} \
			--threads ${task.cpus} \
			--out-dir ${sample_id} \
			--nano-raw $reads"

		//SHELL
    """
    ${flye_command}
    """
}
