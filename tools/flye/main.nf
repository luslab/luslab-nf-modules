#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process flye {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/flye:2.8.1--py38h1c8e9b9_1"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_id}/assembly.fasta"), emit: fasta
        tuple path("${meta.sample_id}/assembly_info.txt"), path("${meta.sample_id}/flye.log"), emit: log

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    flye_command = "flye $args --genome-size ${opts.genome_size} \
			--threads ${task.cpus} \
			--out-dir ${meta.sample_id} \
			--nano-raw $reads"

	//SHELL
    """
    ${flye_command}
    """
}
