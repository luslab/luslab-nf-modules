#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Generate genome index
process star_genome_generate {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    // container 'quay.io/biocontainers/star:2.7.5b--0'
    container 'luslab/nf-modules-star:latest'

    input:
      val opts
      path(fasta)

    output:
      path "genome_index", emit: genome_index

    script:

    // Add the main arguments
    args = "--runMode genomeGenerate --runDirPerm All_RWX --genomeDir genome_index --genomeFastaFiles $fasta "

    // Check and add custom arguments
    if ( opts.args ) {
      if ( opts.args =~ /(--genomeDir)/ ) {
        exit 1, "Error: This module does not support manual setting of --genomeDir. The genome index will appear in ${params.outdir}/star_genome_generate/genome_index. Exit."
      }
      if ( opts.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'genomeGenerate'. You do not need to provide it manually. Exit."
      }
      if ( opts.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the params.modules['star_genome_generate'].args parameter. Exit."
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
        println ("[MODULE] star_genome_generate command: " + index_command)
    }

    // Run read mapping with STAR
    """
    mkdir genome_index
    ${index_command}
    """
}

// Map reads
process star_align_reads {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    // container 'quay.io/biocontainers/star:2.7.5b--0'
    container 'luslab/nf-modules-star:latest'

    input:
      val opts
      tuple val(meta), path(reads)
      path star_index 

    output:
      tuple val(meta), path("*.sam"), optional: true, emit: sam_files
      tuple val(meta), path("*.bam"), path("*.bai"), optional: true, emit: bam_files
      tuple val(meta), path("*.SJ.out.tab"), optional: true, emit: sj_files
      tuple val(meta), path("*.junction"), optional: true, emit: ch_junctions
      tuple val(meta), path("*.ReadsPerGene.out.tab"),  optional: true, emit: reads_per_gene
      tuple val(meta), path("*.Log.final.out"), emit: final_log_files
      tuple val(meta), path("*.Log.out"), emit: out_log_files
      tuple val(meta), path("*.Log.progress.out"), emit: progress_log_files
      path "*.Log.final.out", emit: report

    script:

    // Add the main arguments
    args = "--runMode alignReads --runDirPerm All_RWX --genomeDir $star_index --readFilesIn $reads "

    // Check and add custom arguments
    if ( opts.args ) {
      if ( opts.args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( opts.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( opts.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the params.modules['star_align_reads'].args parameter. Exit."
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

    index_command = "samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam"

    // Construct command line
    map_command = "STAR $args && $index_command"

    // Log
    if (params.verbose) {
        println ("[MODULE] star_align_reads command: " + map_command)
    }

    // Run read mapping with STAR
    """
    ${map_command}
    """
}
