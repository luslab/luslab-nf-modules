#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process seacr {
    publishDir "${params.outdir}/seacr",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                     if (opts.publish_results == "none") null
                     else filename }

    
    // container
    container 'luslab/nf-modules-seacr:latest'

    input:
    val opts
    tuple val(meta), path(bedgraph)
    path control

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "*.log", emit: log

    script:
    outfile_name = seacr_${bedgraph}
    if(opts.outfile_name) {
        outfile_name = opts.outfile_name
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    seacr_command = "SEACR_1.3.sh ${bedgraph} ${control} ${opts.args} ${prefix}"

    """
    ${seacr_command}
    """
}