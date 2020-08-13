#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_intersect_regions {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-bedtools:latest'

    input:
        val opts
        tuple val(meta), path(reads)
        path regions_file

    output:
        tuple val(meta), path("*.bed"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        intersect_command = "bedtools intersect -a ${regions_file} -b $reads ${args} > ${prefix}.bed"
        if (params.verbose){
            println ("[MODULE] bedtools/intersect command: " + intersect_command)
        }

        //SHELL
        """
        ${intersect_command}
        """
}


process bedtools_intersect {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-bedtools:latest'

    input:
        val opts
        tuple val(meta), path(file_a)
        path(file_b)

    output:
        tuple val(meta), path("*.bed"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        intersect_command = "bedtools intersect -a ${file_a} -b ${file_b} ${args} > ${prefix}.bed"
        if (params.verbose){
            println ("[MODULE] bedtools/intersect command: " + intersect_command)
        }

        //SHELL
        """
        ${intersect_command}
        """
}

//Process definition
process bedtools_subtract {
    publishDir "${params.outdir}/bedtools/subtract",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-bedtools:latest'

    input:
        val opts
        tuple val(meta), path(file_a)
        path(file_b)

    output:
        tuple val(meta), path("${prefix}")

    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.bedtools_subtract_args && params.bedtools_subtract_args != '') {
        ext_args = params.bedtools_subtract_args
        args += " " + ext_args.trim()
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    // Construct CL line
    subtract_command = "bedtools subtract -a ${file_a} -b ${file_b} ${args} > '${prefix}'"

    // Log
    if (params.verbose){
        println ("[MODULE] bedtools/subtract command: " + subtract_command)
    }

    //SHELL
    """
    ${subtract_command}
    """
}