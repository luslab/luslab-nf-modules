#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process phyml {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/phyml:3.3.20200621--he7093c6_0"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path ("*_phyml_tree.txt"), emit: tree
        tuple val(meta), path ("*_phyml_stats.txt"), emit: ml_params
        tuple val(meta), path ("*_boot_trees.txt"), optional: true, emit: bootstrap_trees
        tuple val(meta), path ("*_boot_stats.txt"), optional: true, emit: bootstrap_ml_params
        tuple val(meta), path ("*_phyml_trees.txt"), optional: true, emit: multi_trees

    script:

    args = ""
    if(opts.args && opts.args != '') {
        ext_args = opts.args
        args += ext_args.trim()
    }

    phyml_command = "phyml $args --input ${fasta} --datatype ${opts.data_type}"

    if (params.verbose){
        println ("[MODULE] phyml command: " + phyml_command)
    }

    // SHELL
    """
    ${phyml_command}
    """
}
