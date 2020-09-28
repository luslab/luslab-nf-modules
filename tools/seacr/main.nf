#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process seacr {
    publishDir "${params.outdir}/seacr"
    mode: 'copy',
    overwrite: true,
    saveAs: { filename ->
                if (opts.publish_results == "none") null
                else filename }

    
    // container

    input:

    output:

    script:

}