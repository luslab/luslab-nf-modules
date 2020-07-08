#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process guppy_basecalling {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true
    
    input:
        tuple val(sample_id), path(raw)

    output:
        tuple val(sample_id), path("*.fastq"), emit: basecalledSeq

    script:

    // SHELL
    """
    guppy_basecaller
    """
}