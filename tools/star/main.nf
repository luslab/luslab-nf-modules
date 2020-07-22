#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process star_alignReads {
    tag "${sample_id}"

    publishDir "${params.outdir}/star_map",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-star:2.7.5a' // 2.6.1d 2.7.3a

    input:
      tuple val(sample_id), path(reads), path(star_index)

    output:
      tuple val(sample_id), path("*.sam"), optional: true, emit: samFiles
      tuple val(sample_id), path("*.bam"), optional: true, emit: bamFiles
      tuple val(sample_id), path("*.SJ.out.tab"), optional: true, emit: sjFiles
      tuple val(sample_id), path("*.junction"), optional: true, emit: chJunctions
      tuple val(sample_id), path("*.ReadsPerGene.out.tab"),  optional: true, emit: readsPerGene
      tuple val(sample_id), path("*.Log.final.out"), emit: finalLogFiles
      tuple val(sample_id), path("*.Log.out"), emit: outLogFiles
      tuple val(sample_id), path("*.Log.progress.out"), emit: progressLogFiles
      path "*.Log.final.out", emit: report

    script:

    // Initialise argument string
    args = ""

    // Add the main arguments
    args = "--runMode alignReads --genomeDir $star_index --readFilesIn $reads "

    // Check and add custom arguments
    if ( params.star_map_args ) {
      if ( params.star_map_args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( params.star_map_args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( params.star_map_args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the star_map_args parameter. Exit."
      }
      ext_args = params.star_map_args
      args += ext_args.trim() + " "
    }

    // Add the number of threads
    args += "--runThreadN $task.cpus "

    // Add output file name prefix
    args += "--outFileNamePrefix ${sample_id}. "

    // Add compression parameters 
    test_file_name = "$reads"
    if ( "$test_file_name" =~ /(.gz$)/ ) {
      args += "--readFilesCommand gunzip -c "
    } 
    if ( "$test_file_name" =~ /(.bz2$)/ ) {
      args += "--readFilesCommand bunzip2 -c "
    }

    // Add optional input parameters
    if ( params.sjdbGTFfile ) {
      args += "--sjdbGTFfile ${params.sjdbGTFfile} "
    }
    if ( params.sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${params.sjdbFileChrStartEnd} "
    }
    if ( params.varVCFfile ) {
      args += "--varVCFfile ${params.varVCFfile} "
    }

    // Add memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000} " : ''
    avail_mem += task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    args += avail_mem

    // Construct command line
    map_command = "STAR $args"

    // Log
    if (params.verbose) {
        println ("[MODULE] star/map command: " + map_command)
    }

    // Run read mapping with STAR
    """
    ${map_command}
    """
}
