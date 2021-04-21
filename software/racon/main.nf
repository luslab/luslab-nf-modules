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

    container "quay.io/biocontainers/racon:1.4.20--he513fc3_0"

    input:
        val opts
        tuple val(meta), path(reads)
        tuple val(meta), path(overlap_paf)
        tuple val(meta), path(assembly_fasta)

    output:
        tuple val(meta), path("racon_${opts.polish_iter}/${assembly_fasta.simpleName}.fasta"), emit: fasta

    script:
    //Build the command line options
    racon_command = "mkdir -p racon_${opts.polish_iter}; racon --threads ${task.cpus} ${reads} ${overlap_paf} ${assembly_fasta} > racon_${opts.polish_iter}/${assembly_fasta.simpleName}.fasta"

    //SHELL
    """
    ${racon_command}
    """
}
