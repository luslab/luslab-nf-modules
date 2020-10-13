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
    
    // r-base=4.0.2,seacr=1.3,bedtools=2.29.2
    container 'quay.io/biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:5bb5ed4307a8187a7f34730b00431de93688fa59-0'

    input:
    val opts
    tuple val(meta), path(bedgraph)
    tuple val(control_meta), path(control)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    outfile_name = "seacr_${bedgraph}"
    if(opts.outfile_name) {
        outfile_name = opts.outfile_name
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    seacr_command = "SEACR_1.3.sh ${bedgraph} ${control} ${opts.args} ${prefix}"

    if (params.verbose){
        println ("[MODULE] seacr command: " + seacr_command)
    }

    """
    ${seacr_command}
    """
}