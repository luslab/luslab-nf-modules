#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process infernal_cmscan {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/infernal:1.1.3--h516909a_0"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.cmscan"), emit: cmscan
        tuple val(meta), path("*.tbl"), emit: tbl

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    cmscan_command = "cmscan $args --cpu ${task.cpus} --tblout ${fasta.simpleName}.tbl -Z ${opts.search_space} --clanin ${opts.clanin} ${opts.db} ${fasta} > ${fasta.simpleName}.cmscan"

    if (params.verbose){
        println ("[MODULE] cmscan command: " + cmscan_command)
    }
    //SHELL
    """
    ${cmscan_command}
    """
}
