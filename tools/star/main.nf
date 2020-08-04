#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Generate genome index
process star_genomeGenerate {
    publishDir "${params.outdir}/star_genomeGenerate",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-star:2.7.5a'

    input:
      path(fasta)

    output:
      path "genome_index", emit: genomeIndex
      path "genome_index/chrName.txt", emit: chrNameFile
      path "genome_index/Log.out", emit: report

    script:

    // Initialise argument string
    args = ""

    // Add the main arguments
    args = "--runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles $fasta "

    // Check and add custom arguments
    if ( params.star_genomeGenerate_args ) {
      if ( params.star_genomeGenerate_args =~ /(--genomeDir)/ ) {
        exit 1, "Error: This module does not support manual setting of --genomeDir. The genome index will appear in ${params.outdir}/star_genomeGenerate/genome_index. Exit."
      }
      if ( params.star_genomeGenerate_args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'genomeGenerate'. You do not need to provide it manually. Exit."
      }
      if ( params.star_genomeGenerate_args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the star_genomeGenerate_args parameter. Exit."
      }
      ext_args = params.star_genomeGenerate_args
      args += ext_args.trim() + " "
    }

    // Add the number of threads
    args += "--runThreadN $task.cpus "

    // Add optional input parameters
    if ( params.star_genomeGenerate_sjdbGTFfile ) {
      args += "--sjdbGTFfile ${params.star_genomeGenerate_sjdbGTFfile} "
    }
    if ( params.star_genomeGenerate_sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${params.star_genomeGenerate_sjdbFileChrStartEnd} "
    }

    // Add memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000} " : ''
    args += avail_mem

    // Construct command line
    index_command = "STAR $args"

    // Log
    if (params.verbose) {
        println ("[MODULE] star_genomeGenerate command: " + index_command)
    }

    // Run read mapping with STAR
    """
    mkdir genome_index
    ${index_command}
    """
}

// Map reads
process star_alignReads {
    tag "${sample_id}"

    publishDir "${params.outdir}/star_alignReads",
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
    if ( params.star_alignReads_args ) {
      if ( params.star_alignReads_args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( params.star_alignReads_args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( params.star_alignReads_args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the star_alignReads_args parameter. Exit."
      }
      ext_args = params.star_alignReads_args
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
    if ( params.star_alignReads_sjdbGTFfile ) {
      args += "--sjdbGTFfile ${params.star_alignReads_sjdbGTFfile} "
    }
    if ( params.star_alignReads_sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${params.star_alignReads_sjdbFileChrStartEnd} "
    }
    if ( params.star_alignReads_varVCFfile ) {
      args += "--varVCFfile ${params.star_alignReads_varVCFfile} "
    }

    // Add memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000} " : ''
    avail_mem += task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    args += avail_mem

    // Construct command line
    map_command = "STAR $args"

    // Log
    if (params.verbose) {
        println ("[MODULE] star_alignReads command: " + map_command)
    }

    // Run read mapping with STAR
    """
    ${map_command}
    """
}
