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
    container 'luslab/nf-modules-seacr:latest'

    input:
    val opts
    tuple val(meta), path(bedgraph)
    path(control)

    output:
    tuple val(meta), path("*.bed")

    script:

    seacr_command = "SEACR_1.3.sh ${bedgraph} ${control} ${opts.args} ${opts.out_file}.bed"

}