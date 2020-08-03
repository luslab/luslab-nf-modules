#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_intersect {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-hisat2:latest'

    input:
        val opts
        tuple val(meta), path(reads)
        path spice_sites
        path inexed_genome

    output:
        tuple val(meta), path("*.sam"), emit: sam

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        //SHELL
        readList = reads.collect{it.toString()}
        if (readList.size > 1){
            hisat2_command = "hisat2 -x ${inexed_genome} -p ${task.cpus} --known-splicesite-infile ${spice_sites} \
            --met-file ${prefix}.txt -U ${reads} -S ${prefix}.sam ${opts.args}"
        } else {
            hisat2_command = "hisat2 -x ${inexed_genome} -p ${task.cpus} --known-splicesite-infile ${spice_sites} \
            --met-file ${prefix}.txt -1 ${reads[0]} -2 ${reads[1]} -S ${prefix}.sam ${opts.args}"
        }

        if (params.verbose){
            println ("[MODULE] hisat2 command: " + hisat2_command)
        }

        //SHELL
        """
        ${hisat2_command}
        """
}