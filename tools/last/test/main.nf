#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for last...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {last} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/test-nanopore.fastq.gz"],
]

// Define regions file input channel
Channel
    .fromPath("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_test_fasta}

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_fastq}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run last
    bedtools_intersect (params.modules['last'], ch_test_fastq, ch_test_fasta )

    // Collect file names and view output
    last.out.mappedReads | view
}