#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Samtools index
process samtools_index {
    publishDir "${params.outdir}/samtools/index",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*.bam.bai"), emit: baiFiles
 
    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.samtools_index_args && params.samtools_index_args != '') {
        ext_args = params.samtools_index_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    index_command = "samtools index ${args} -@ ${task.cpus} ${bam[0]}"

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
    publishDir "${params.outdir}/samtools/view",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        tuple val(sample_id), path(bam_bai)

    output:
        tuple val(sample_id), path("*.b*"), emit: bamFiles
 
    script:

    // Check main args string exists and strip whitespace
    args = "-b -h"
    if(params.samtools_view_args && params.samtools_view_args != '') {
        ext_args = params.samtools_view_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    view_command = "samtools view ${args} -@ ${task.cpus} -o ${bam_bai[0].simpleName}.filt.bam ${bam_bai[0]} ${params.samtools_view_region}"

    // Log
    if (params.verbose){
        println ("[MODULE] samtools/view command: " + view_command)
    }
    
    """
    ${view_command}
    samtools index -@ ${task.cpus} ${bam_bai[0].simpleName}.filt.bam
    """
}

// Samtools faidx indexes fasta files
process samtools_faidx {
    publishDir "${params.outdir}/samtools/faidx",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        path fasta

    output:
        tuple path("*.fa", includeInputs: true), path("*.fai"), emit: indexedFiles
 
    script:

    // Construct CL line
    command = "samtools faidx $fasta"

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
        tuple val(meta), path("*.bam"), emit: bam

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        sort_command = "samtools sort -0 $reads ${args} > ${prefix}.bed"
        if (params.verbose){
            println ("[MODULE] samtools/sort command: " + sort_command)
        }

        //SHELL
        """
        ${sort_command}
        """
}