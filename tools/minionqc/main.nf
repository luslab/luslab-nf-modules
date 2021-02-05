#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process minionqc {
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-minionqc:base-1.0.2"

    input:
        val opts
        tuple val(meta), path(sequencing_summary)

    output:
        tuple val(meta), path("minionqc"), emit: minionqc_output_dir

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    minionqc_command = "MinIONQC.R -p ${task.cpus} ${args} -o minionqc -i $sequencing_summary"

    if (params.verbose){
        println ("[MODULE] minionqc command: " + minionqc_command)
    }

    //SHELL
    """
    ${minionqc_command}
    """
}
