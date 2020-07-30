#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for minimap2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {minimap2} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testData = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/test-nanopore.fastq.gz"],
]

//Define test data input channel
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
    // Run minimap2
    minimap2 (params.modules['minimap2'], ch_test_fastq, ch_test_fasta )

    // Collect and view output
    minimap2.out.mappedReads | view
}