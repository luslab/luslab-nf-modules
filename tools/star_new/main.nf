#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// --outfile name prefix
// Set to sample ID
params.outfile_prefix_sampleid = true

// Switch for paired-end reads 
params.paired_end = false

// Process definition
process star_map {
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
    args = ""

    // Set the main arguments
    if (params.paired_end) {
      args = "--genomeDir $star_index --readFilesIn ${reads[0]} ${reads[1]} "
    } else {
      args = "--genomeDir $star_index --readFilesIn $reads "
    }

    // Combining the custom arguments and creating star args
    if (params.star_map_args) {
      ext_args = params.star_map_args
      args += " " + ext_args.trim()
    }

    // Set the number of threads
    args += "--runThreadN $task.cpus "

    // Output file name prefix
    if (params.outfile_prefix_sampleid) {
      args += "--outFileNamePrefix ${sample_id}. "
    }

    // Compression parameters 
    if (params.paired_end) {
      test_file_name = "${reads[0]}"
    } else {
      test_file_name = "$reads"
    }

    if ("$test_file_name" =~ /(.gz$)/) {
      args += "--readFilesCommand gunzip -c "
    } 
    if ("$test_file_name" =~ /(.bz2$)/) {
      args += "--readFilesCommand bunzip2 -c "
    }

    // Set memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000} " : ''
    avail_mem += task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    args += avail_mem

    // Construct command line
    map_command = "STAR $arg"

    // Log
    if (params.verbose) {
        println ("[MODULE] star/map command: " + map_command)
    }

    // Run read mapping with STAR
    """
    ${map_command}
    """
}
