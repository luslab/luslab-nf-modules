#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "star/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'map'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'map'

// STAR parameters
params.internal_custom_args = ''

//--outfile name prefix
// Set to sample ID
params.internal_outfile_prefix_sampleid = true

//Switch for paired-end files 
params.internal_paired_end = false

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

// Trimming reusable component
process map {
    tag "${sample_id}"

    publishDir "star/map/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

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

    shell:
    
    // Set the main arguments
    if (params.internal_paired_end){
      star_args = "--genomeDir $star_index --readFilesIn ${reads[0]} ${reads[1]} "
    } else {
      star_args = "--genomeDir $star_index --readFilesIn $reads "
    }

    // Combining the custom arguments and creating star args
    if (params.internal_custom_args){
      star_args += "$params.internal_custom_args "
    }

    //RunThread param
     star_args += "--runThreadN $task.cpus "

    //outfile name prefix
    if (params.internal_outfile_prefix_sampleid){
      star_args += "--outFileNamePrefix ${sample_id}. "
    }

    // Compression parameters 
    if (params.internal_paired_end) {
      test_file_name = "${reads[0]}"
    } else {
      test_file_name = "$reads"
    }

    if ("$test_file_name" =~ /(.gz$)/){
      star_args += "--readFilesCommand gunzip -c "
    } 
    if ("$test_file_name" =~ /(.bz2$)/){
      star_args += "--readFilesCommand bunzip2 -c "
    }

    // Set memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    avail_mem += task.memory ? " --limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    star_args += avail_mem
    
    """
    STAR $star_args
    """
}