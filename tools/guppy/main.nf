#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process guppy_basecaller {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true
    
    container "luslab/nf-modules-guppy:latest"

    input:
        path(reads)

    output:
        tuple path("*.fastq"), path("*.log"), path("*.txt"), path("*.js"), emit: basecalledSeq
        

    script:

    // SHELL
    """
    guppy_basecaller --input_path $reads --save_path . --flowcell ${params.guppy_flowcell} --kit ${params.guppy_kit}
    """
}