#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_intersect_regions {
    tag "${meta.sample_id}" 

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

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

        intersect_command = "bedtools intersect ${args} -a ${regions_file} -b $reads > ${prefix}"
        if (params.verbose){
            println ("[MODULE] bedtools/intersect_regions command: " + intersect_command)
        }

        //SHELL
        """
        ${intersect_command}
        """
}


process bedtools_intersect {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(file_a)
        path(file_b)

    output:
        tuple val(meta), path("${prefix}"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        intersect_command = "bedtools intersect -a ${file_a} -b ${file_b} ${args} > ${prefix}"
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
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(file_a)
        path(file_b)

    output:
        tuple val(meta), path("${prefix}"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        subtract_command = "bedtools subtract -a ${file_a} -b ${file_b} ${args} > ${prefix}"
        if (params.verbose){
            println ("[MODULE] bedtools/subtract command: " + subtract_command)
        }

        //SHELL
        """
        ${subtract_command}
        """
}

//Process definition
process bedtools_bamtobed {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(bam_file), path(bai_file)

    output:
        tuple val(meta), path("${prefix}"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        bamtobed_command = "bedtools bamtobed -i ${bam_file} ${args} > ${prefix}"
        if (params.verbose){
            println ("[MODULE] bedtools/bamtobed command: " + bamtobed_command)
        }

    //SHELL
    """
    ${bamtobed_command}
    """

}

//Process definition
process bedtools_genomecov {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(aligned_file)
        path genome 

    output:
        tuple val(meta), path("${prefix}"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        genomecov_command = "bedtools genomecov ${args} -i ${aligned_file} -g ${genome} > ${prefix}"
        if (params.verbose){
            println ("[MODULE] bedtools/genomecov command: " + genomecov_command)
        }

    //SHELL
    """
    ${genomecov_command}
    """

}

process bedtools_genomecov_scale {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(aligned_file)
        path genome
        val scale

    output:
        tuple val(meta), path("${prefix}"), emit: bedgraph

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        if (args == ''){
            genomecov_command = "bedtools genomecov -i ${aligned_file} -g ${genome} -bg -scale ${scale} > ${prefix}"
        } else {
            genomecov_command = "bedtools genomecov ${args} -i ${aligned_file} -g ${genome} -bg -scale ${scale} > ${prefix}"
        }
        if (params.verbose){
            println ("[MODULE] bedtools/genomecov command: " + genomecov_command)
        }

    //SHELL
    """
    ${genomecov_command}
    """

}

process bedtools_genomecov_bam {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(bam_file), path(bai_file)

    output:
        tuple val(meta), path("${prefix}"), emit: bed

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        genomecov_command = "bedtools genomecov ${args} -ibam ${bam_file} > ${prefix}"
        if (params.verbose){
            println ("[MODULE] bedtools/genomecov command: " + genomecov_command)
        }

    //SHELL
    """
    ${genomecov_command}
    """

}

process bedtools_genomecov_scale_bam {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        val opts
        tuple val(meta), path(bam), path(bai)
        val scale

    output:
        tuple val(meta), path("${prefix}"), emit: bedgraph

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        if (args == ''){
            genomecov_command = "bedtools genomecov -ibam ${bam} -bg -scale ${scale} > ${prefix}"
        } else {
            genomecov_command = "bedtools genomecov ${args} -ibam ${bam} -bg -scale ${scale} > ${prefix}"
        }
        if (params.verbose){
            println ("[MODULE] bedtools/genomecov command: " + genomecov_command)
        }

    //SHELL
    """
    ${genomecov_command}
    """

}