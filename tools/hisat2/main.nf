#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hisat2_build {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${genome}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2'

    input:
        val opts
        path genome

    output:
        path "${opts.index_dir}", emit: genome_index

    script:

        //SHELL
        """
        mkdir ${opts.index_dir}
        hisat2-build -p ${task.cpus} ${genome} ${opts.index_dir}/${opts.prefix}
        """
}

process hisat2_splice_sites {
    label "min_cores"
    label "min_mem"
    label "regular_queue"

    tag "${gtf}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2'

    input:
        val opts
        path gtf

    output:
        path "${opts.prefix}${opts.suffix}", emit: splice_sites

    script:

        //SHELL
        """
        mkdir hisat2_index
        hisat2_extract_splice_sites.py ${gtf} > ${opts.prefix}${opts.suffix}
        """
}
 

process hisat2_splice_align {
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

    container 'quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2'

    input:
        val opts
        tuple val(meta), path(reads)
        path genome_index
        path splice_sites

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
            hisat2_align_command = "hisat2 -x ${genome_index}/index -p ${task.cpus} --known-splicesite-infile ${splice_sites} \
            --met-file ${prefix}.txt -1 ${reads[0]} -2 ${reads[1]} -S ${prefix}.sam ${opts.args}"
        } else {
            hisat2_align_command = "hisat2 -x ${genome_index}/index -p ${task.cpus} --known-splicesite-infile ${splice_sites} \
            --met-file ${prefix}.txt -U ${reads} -S ${prefix}.sam ${opts.args}"
        }

        if (params.verbose){
            println ("[MODULE] hisat2 command: " + hisat2_align_command)
        }

        //SHELL
        """
        ${hisat2_align_command}
        """
}

process hisat2_align {
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

    container 'quay.io/biocontainers/hisat2:2.2.1--he1b5a44_2'

    input:
        val opts
        tuple val(meta), path(reads)
        path genome_index

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
            hisat2_align_command = "hisat2 -x ${genome_index}/index -p ${task.cpus} \
            --met-file ${prefix}.txt -1 ${reads[0]} -2 ${reads[1]} -S ${prefix}.sam ${opts.args}"
        } else {
            hisat2_align_command = "hisat2 -x ${genome_index}/index -p ${task.cpus} \
            --met-file ${prefix}.txt -U ${reads} -S ${prefix}.sam ${opts.args}"
        }

        if (params.verbose){
            println ("[MODULE] hisat2 command: " + hisat2_align_command)
        }

        //SHELL
        """
        ${hisat2_align_command}
        """
}