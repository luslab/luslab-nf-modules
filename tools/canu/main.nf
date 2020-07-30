#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process minionqc {
    publishDir "${params.outdir}/minionqc",
        mode: "copy", overwrite: true

    container "luslab/nf-modules-minionqc"

    input:
        tuple val(sample_id), path(sequencing_summary)

    output:
        tuple val(sample_id), path("*_minionqc"), emit: minionqcOutputs

    script:

    minionqc_command = "Rscript /MinIONQC.R -p ${task.cpus} -o ${sample_id}_minionqc -f pdf -i $sequencing_summary"

		//SHELL
    """
    ${minionqc_command}
    """
}
