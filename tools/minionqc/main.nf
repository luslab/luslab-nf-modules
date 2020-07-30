#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process minionqc {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-minionqc:latest"

    input:
        val opts
        tuple val(meta), path(sequencing_summary)

    output:
        tuple val(meta), path("*_minionqc"), emit: minionqcOutputs

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    minionqc_command = "Rscript MinIONQC.R -p ${task.cpus} -o ${prefix}_minionqc -i $sequencing_summary"

    if (params.verbose){
        println ("[MODULE] minionqc command: " + minionqc_command)
    }

	//SHELL
    """
    ${minionqc_command}
    """
}