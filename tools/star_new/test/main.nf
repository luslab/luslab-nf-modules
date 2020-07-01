#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for STAR mapping...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.star_map_args = '--outFilterMultimapNmax 20'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include star_map as map_se from '../main.nf'
include star_map as map_pe from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Single-end test reads
testMetaDataSingleEnd = [
  ['Sample1', "$baseDir/input/prpf8_eif4a3_rep1.Unmapped.fq"],
  ['Sample2', "$baseDir/input/prpf8_eif4a3_rep2.Unmapped.fq"]
]

// Paired-end test reads
testMetaDataPairedEnd = [
  ['Sample1', "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"]
]

// Channel for single-end reads 
 Channel
    .from( testMetaDataSingleEnd )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData_single_end }

// Channel for paired-end reads
  Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData_paired_end }

// Run tests
workflow {
    // Run dedup
    map_se ( ch_testData_single_end )

    // Collect file names and view output for single-end read mapping
    ...

    // Collect file names and view output for single-end read mapping
    ...
}