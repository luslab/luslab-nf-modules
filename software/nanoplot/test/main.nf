#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

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
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_nanopore_reads = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_all.fastq.gz"],
]

//Define test data input channels
Channel
    .from(test_nanopore_reads)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_nanoporeread_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    nanoplot(params.modules["nanoplot"], ch_nanoporeread_data)

    // Collect file names and view output
    nanoplot.out.report | view

    assert_channel_count( nanoplot.out.report, "output", 1)
}
