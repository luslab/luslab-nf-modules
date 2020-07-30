#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process minimap2 {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true

    container "luslab/nf-modules-minimap2:latest"

    input:

    output:

    script:
    """
    minimap2 
    """
}