#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for porechop...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {porechop} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_nanopore = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_all.fastq.gz"],
]

//Define test data input channels

//Single end
Channel
    .from(test_nanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    porechop ( params.modules["porechop"], ch_fastq )

    // Collect file names and view output
    porechop.out.fastq | view

    // Verify channel counts
    assert_channel_count(porechop.out.fastq, "porechop'd fastq", 1)
}
