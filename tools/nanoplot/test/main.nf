#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Log out
log.info ("Starting tests for NanoPlot...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {nanoplot} from "../main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testNanoporeReads = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/fastq/test-nanopore.fastq.gz"],
]

//Define test data input channels
Channel
    .from(testNanoporeReads)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_nanoporeread_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    nanoplot(params.modules['nanoplot'], ch_nanoporeread_data)

    // Collect file names and view output
    nanoplot.out.nanoplotOutputs | view
}
