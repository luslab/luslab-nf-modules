#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process umitools_dedup {
    tag "${meta.sample_id}"
    
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    // Biocontainer with umi_tools, samtools
    container 'quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:b4c5bc18f99e35cd061328b98f501aa7b9717308-0'

    input:
        val(opts)
        tuple val(meta), path(bam), path(bai)
       
    output:
        tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: dedupBam
        path "*.log", emit: report

    script:

    // Init
    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    args = "--log=${prefix}.log "

    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    // Construct CL line
    dedup_command = "umi_tools dedup ${args} -I ${bam[0]} -S ${prefix}.bam --output-stats=${prefix}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    samtools index ${prefix}.bam
    """
}
