#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process raven {
    label "high_cores"
    label "max_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/raven-assembler:1.4.0--h8b12597_0"

    input:
        val opts
        // Note that "reads" can be FASTA or FASTQ format and compressed or not.
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_id}/assembly.fasta"), emit: fasta
        tuple val(meta), path("${meta.sample_id}/assembly_graph.gfa"), emit: gfa

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    raven_command = "mkdir -p ${meta.sample_id} && raven $args --threads ${task.cpus} --polishing-rounds ${opts.polishing_rounds} --graphical-fragment-assembly ${meta.sample_id}/assembly_graph.gfa ${reads} > ${meta.sample_id}/assembly.fasta"

    //SHELL
    """
    ${raven_command}
    """
}
