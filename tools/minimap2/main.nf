#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process minimap2_index {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"

    input:
        val opts
        tuple val(meta), path (fasta)

    output:
        tuple val(meta), path ("*.mmi"), emit: mmi

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    minimap2_command = "minimap2 $args -t ${task.cpus} -d ${fasta.simpleName}.mmi ${fasta}"

    if (params.verbose){
        println ("[MODULE] minimap2_index command: " + minimap2_command)
    }

    // SHELL
    """
    ${minimap2_command}
    """
}

process minimap2_paf {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"

    input:
        val opts
        tuple val(meta), path (fasta)
        tuple val(meta), path (reads)


    output:
        tuple val(meta), path ("*.paf"), emit: paf

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    minimap2_command = "minimap2 $args -t ${task.cpus} ${fasta} ${reads} > ${fasta.simpleName}.paf"

    if (params.verbose){
        println ("[MODULE] minimap2_paf command: " + minimap2_command)
    }

    // SHELL
    """
    ${minimap2_command}
    """
}


process minimap2_sam {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"

    input:
        val opts
        tuple val(meta), path (fasta)
        tuple val(meta), path (reads)

    output:
        tuple val(meta), path ("*.sam"), emit: sam

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    minimap2_command = "minimap2 $args -t ${task.cpus} -a ${fasta} ${reads} > ${fasta.simpleName}.sam"

    if (params.verbose){
        println ("[MODULE] minimap2_sam command: " + minimap2_command)
    }

    // SHELL
    """
    ${minimap2_command}
    """
}
