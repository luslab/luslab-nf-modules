#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for samtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.samtools_index_args = ''
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include samtools_index from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    ['Sample1', "$baseDir/input/sample1.bam"],
    ['Sample2', "$baseDir/input/sample2.bam"]
]

Channel
    .from( testData )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testData }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run star
    samtools_index ( ch_testData )

    // Collect file names and view output
    samtools_index.out.baiFiles | view
}
