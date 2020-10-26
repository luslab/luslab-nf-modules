#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Samtools index
process samtools_index {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_id}.bam", includeInputs: true), path("${prefix}"), emit: bam
 
    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}.bai" : "${meta.sample_id}.bai"

        index_command = "samtools index ${args} -@ ${task.cpus} ${reads} > ${prefix}"

        if (params.verbose){
            println ("[MODULE] samtools/index command: " + index_command)
        }
        
    """
    ${index_command}
    mv ${reads} ${meta.sample_id}.bam
    """
}

// Samtools view - only works for bam files - requires bam and bai
process samtools_view {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

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
    view_command = "samtools view ${args} -@ ${task.cpus} ${reads} ${samtools_view_region} > ${prefix}"

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
    tag "${fasta}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

    input:
        val opts
        path fasta

    output:
        tuple path("*.fa", includeInputs: true), path("*.fai"), emit: fasta
        path "*.fai", emit: fai
 
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
    tag "${meta.sample_id}"
    
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'
    
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

        sort_command = "samtools sort ${args} -o ${prefix} $reads"
        if (params.verbose){
            println ("[MODULE] samtools/sort command: " + sort_command)
        }

        //SHELL
        """
        ${sort_command}
        """
}
