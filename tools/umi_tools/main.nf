#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process umitools_dedup {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-umi_tools:latest'

    input:
        val(opts)
        tuple val(meta), path(bam)
       
    output:
        tuple val(meta), path("${prefix}${opts.suffix}.bam"), path("${prefix}${opts.suffix}.bam.bai"), emit: dedupBam
        path "*.log", emit: report

    script:

    // Init
    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    args = "--log=${prefix}${opts.suffix}.log "

    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    // Construct CL line
    dedup_command = "umi_tools dedup ${args} -I ${bam[0]} -S ${prefix}${opts.suffix}.bam --output-stats=${prefix}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    samtools index ${prefix}${opts.suffix}.bam
    """
}
