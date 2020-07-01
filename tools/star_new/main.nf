#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// --outfile name prefix
// Set to sample ID
params.outfile_prefix_sampleid = true

// Switch for paired-end reads 
params.paired_end = false

// Process definition
process map {
    tag "${sample_id}"

    publishDir "${params.outdir}/star/map",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-star:latest'

    input:
      tuple val(sample_id), path(reads), path(star_index)

    output:
      tuple val(sample_id), path("*.sam"), optional: true, emit: samFiles
      tuple val(sample_id), path("*.bam"), optional: true, emit: bamFiles
      tuple val(sample_id), path("*SJ.out.tab"), emit: sjFiles
      tuple val(sample_id), path("*Log.final.out"), emit: finalLogFiles
      tuple val(sample_id), path("*Log.out"), emit: outLogFiles
      tuple val(sample_id), path("*Log.progress.out"), emit: progressLogFiles
      path "*Log.final.out", emit: report

    script:

    // Init
    args = "--log=${sample_id}.dedup.log"

    // Check main args string exists and strip whitespace
    if(params.umitools_dedup_args) {
        ext_args = params.umitools_dedup_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    dedup_command = "umi_tools dedup ${args} -I ${bam[0]} -S ${sample_id}.dedup.bam --output-stats=${sample_id}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    samtools index ${sample_id}.dedup.bam
    """
}
