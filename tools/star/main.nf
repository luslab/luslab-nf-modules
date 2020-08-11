#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Generate genome index
process star_genomeGenerate {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-star:2.7.5a'

    input:
      val opts
      path(fasta)

    output:
      path "genome_index", emit: genomeIndex
      path "genome_index/chrName.txt", emit: chrNameFile
      path "genome_index/Log.out", emit: report

    script:

    // Add the main arguments
    args = "--runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles $fasta "

    // Check and add custom arguments
    if ( opts.args ) {
      if ( opts.args =~ /(--genomeDir)/ ) {
        exit 1, "Error: This module does not support manual setting of --genomeDir. The genome index will appear in ${params.outdir}/star_genomeGenerate/genome_index. Exit."
      }
      if ( opts.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'genomeGenerate'. You do not need to provide it manually. Exit."
      }
      if ( opts.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the star_genomeGenerate_args parameter. Exit."
      }
      ext_args = opts.args
      args += ext_args.trim() + " "
    }

    // Add the number of threads
    args += "--runThreadN $task.cpus "

    // Add optional input parameters
    if ( opts.sjdbGTFfile ) {
      args += "--sjdbGTFfile ${opts.sjdbGTFfile} "
    }
    if ( opts.sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${opts.sjdbFileChrStartEnd} "
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
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-star:2.7.5a' // 2.6.1d 2.7.3a

    input:
      val opts
      tuple val(meta), path(reads)
      path star_index 

    output:
      tuple val(meta), path("*.sam"), optional: true, emit: samFiles
      tuple val(meta), path("*.bam"), optional: true, emit: bamFiles
      tuple val(meta), path("*.SJ.out.tab"), optional: true, emit: sjFiles
      tuple val(meta), path("*.junction"), optional: true, emit: chJunctions
      tuple val(meta), path("*.ReadsPerGene.out.tab"),  optional: true, emit: readsPerGene
      tuple val(meta), path("*.Log.final.out"), emit: finalLogFiles
      tuple val(meta), path("*.Log.out"), emit: outLogFiles
      tuple val(meta), path("*.Log.progress.out"), emit: progressLogFiles
      path "*.Log.final.out", emit: report

    script:

    // Add the main arguments
    args = "--runMode alignReads --genomeDir $star_index --readFilesIn $reads "

    // Check and add custom arguments
    if ( opts.args ) {
      if ( opts.args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( opts.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( opts.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the star_alignReads_args parameter. Exit."
      }
      ext_args = opts.args
      args += ext_args.trim() + " "
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    // Add the number of threads
    args += "--runThreadN $task.cpus "

    // Add output file name prefix
    args += "--outFileNamePrefix ${prefix}. "

    // Add compression parameters 
    test_file_name = "$reads"
    if ( "$test_file_name" =~ /(.gz$)/ ) {
      args += "--readFilesCommand gunzip -c "
    } 
    if ( "$test_file_name" =~ /(.bz2$)/ ) {
      args += "--readFilesCommand bunzip2 -c "
    }

    // Add optional input parameters
    if ( opts.sjdbGTFfile ) {
      args += "--sjdbGTFfile ${opts.sjdbGTFfile} "
    }
    if ( opts.sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${opts.sjdbFileChrStartEnd} "
    }
    if ( opts.varVCFfile ) {
      args += "--varVCFfile ${opts.varVCFfile} "
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
