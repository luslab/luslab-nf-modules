#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl = 2

// Process definition
process busco_genome {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "ezlabgva/busco:v4.1.4_cv1"

    input:
        val opts
        tuple val(meta), path(genome)

    output:
        tuple val(meta), path("busco"), emit: report

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    busco_command = "AUGUSTUS_CONFIG_PATH=augustus_config busco $args -m genome -c ${task.cpus} -i $genome -o busco"

    if (params.verbose){
        println ("[MODULE] busco command: " + busco_command)
    }

	//SHELL
    """
    cp -a /augustus/config  augustus_config
    ${busco_command}
    """
}
