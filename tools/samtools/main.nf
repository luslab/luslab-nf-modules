#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Samtools index
process samtools_index {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-samtools:latest'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${prefix}"), emit: bai
 
    script:

        // Check main args string exists and strip whitespace
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        // Construct CL line
        index_command = "samtools index ${args} -@ ${task.cpus} ${reads[0]}"

        // Log
        if (params.verbose){
            println ("[MODULE] samtools/index command: " + index_command)
        }
        
    """
    ${index_command}
    """
}


// Samtools view - only works for bam files - requires bam and bai
process samtools_view {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-samtools:latest'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(prefix), emit: bam
 
    script:
        // Check main args string exists and strip whitespace
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        samtools_view_region = ""
        if(opts.samtools_view_region && opts.samtools_view_region != '') {
            samtools_view_region = opts.samtools_view_region
        }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"


    // Construct CL line
    view_command = "samtools view ${args} -@ ${task.cpus} -o ${prefix} ${reads} ${samtools_view_region}"

    // Log
    if (params.verbose){
        println ("[MODULE] samtools/view command: " + view_command)
    }
    
    """
    ${view_command}
    """
}

// Samtools faidx indexes fasta files
process samtools_faidx {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-samtools:latest'

    input:
        val opts
        path fasta

    output:
        tuple path("*.fa", includeInputs: true), path("*.fai"), emit: indexedFasta
 
    script:
        // Check main args string exists and strip whitespace
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        // Construct CL line
        command = "samtools faidx ${args} $fasta"

        // Log
        if (params.verbose){
            println ("[MODULE] samtools/faidx command: " + command)
        }
        
    """
    ${command}
    """
}


process samtools_sort {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-samtools:latest'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${prefix}"), emit: bam

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        sort_command = "samtools sort ${args} -o $reads > ${prefix}"
        if (params.verbose){
            println ("[MODULE] samtools/sort command: " + sort_command)
        }

        //SHELL
        """
        ${sort_command}
        """
}


// Samtools view
process samtools_view2 {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-samtools:latest'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(prefix), emit: bam
 
    script:

        // Check main args string exists and strip whitespace
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        // Construct CL line
        view_command = "samtools view ${args} -@ ${task.cpus}"

        // Log
        if (params.verbose){
            println ("[MODULE] samtools/view command: " + view_command)
        }
        
    """
    ${view_command}
    """
}