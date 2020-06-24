#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting STAR mapping module")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

/*include map as paired_end_map from '../main.nf' addParams(map_custom_args: 
      "--genomeLoad NoSharedMemory \
      --outFilterMultimapNmax 1 \
      --outFilterMultimapScoreRange 1 \
      --outSAMattributes All \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterType BySJout \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --outFilterScoreMin 10  \
      --alignEndsType Extend5pOfRead1 \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate", map_paired_end: true, map_process_name: 'paired_end_map')*/

include map as single_end_map from '../main.nf' addParams(map_custom_args: 
      "--genomeLoad NoSharedMemory \
      --outFilterMultimapNmax 1 \
      --outFilterMultimapScoreRange 1 \
      --outSAMattributes All \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterType BySJout \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --outFilterScoreMin 10  \
      --alignEndsType Extend5pOfRead1 \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate", map_process_name: 'single_end_map')

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

params.genome_index = "$baseDir/input/reduced_star_index"


//test data for single-end reads
testMetaData = [
  ['Sample1', "$baseDir/input/prpf8_eif4a3_rep1.Unmapped.fq"],
  ['Sample2', "$baseDir/input/prpf8_eif4a3_rep2.Unmapped.fq"]
]

//test data for paired-end reads
testMetaDataPairedEnd = [
  ['Sample1', "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"],
]

//Channel for single-end reads 
 Channel
    .from( testMetaData )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData }

// Channel for paired-end reads
  Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData_paired_end }

/*------------------------------------------------------------------------------------*/

// Run workflow
// Choose between single-end and paired-end reads
workflow {
    // Run map on paired-end reads
    /*paired_end_map( ch_testData_paired_end )

    // Collect file names and view output
    paired_end_map.out.bamFiles.collect() | view // // custom args changes default output to bam files instead of sam files
    //paired_end_map.out.samFiles.collect() | view 
    paired_end_map.out.sjFiles.collect() | view
    paired_end_map.out.finalLogFiles.collect() | view
    paired_end_map.out.outLogFiles.collect() | view
    paired_end_map.out.progressLogFiles.collect() | view */
    
    // Run map on single-end reads
    single_end_map( ch_testData )

    // Collect file names and view output
    single_end_map.out.bamFiles.collect() | view // custom args changes default output to bam files instead of sam files
    //single_end_map.out.samFiles.collect() | view -> default output
    single_end_map.out.sjFiles.collect() | view
    single_end_map.out.finalLogFiles.collect() | view
    single_end_map.out.outLogFiles.collect() | view
    single_end_map.out.progressLogFiles.collect() | view 
}