#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for Shasta...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {shasta} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_nanopore = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq"],
]

//Define test data input channels

//Single end
Channel
    .from(test_data_nanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    shasta ( params.modules["shasta"], ch_fastq_data )

    // Collect file names and view output
    shasta.out.fasta | view
    shasta.out.gfa | view
    shasta.out.report | view

    // Verify channel counts
    assert_channel_count(shasta.out.fasta, "assembly FASTA", 1)
    assert_channel_count(shasta.out.gfa, "assembly graph", 1)
    assert_channel_count(shasta.out.report, "log files", 1)
}
