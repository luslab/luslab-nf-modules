#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process cdhit_prot {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/cd-hit:4.8.1--h8b12597_3"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.fasta"), emit: fasta
        tuple val(meta), path("*.clstr"), emit: clstr

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    cdhit_command = "cd-hit $args -T ${task.cpus} -c ${opts.identity} -i ${fasta} -o ${fasta.simpleName}-${opts.identity}.fasta -p ${fasta.simpleName}-${opts.identity}.clstr"

    if (params.verbose){
        println ("[MODULE] cdhit_prot command: " + cdhit_command)
    }

	//SHELL
    """
    ${cdhit_command}
    """
}

process cdhit_nucl {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/cd-hit:4.8.1--h8b12597_3"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.fasta"), emit: fasta
        tuple val(meta), path("*.clstr"), emit: clstr

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    cdhit_command = "cd-hit-est $args -T ${task.cpus} -c ${opts.identity} -i ${fasta} -o ${fasta.simpleName}-${opts.identity}.fasta -p ${fasta.simpleName}-${opts.identity}.clstr"

    if (params.verbose){
        println ("[MODULE] cdhit_nucl command: " + cdhit_command)
    }

	//SHELL
    """
    ${cdhit_command}
    """
}
