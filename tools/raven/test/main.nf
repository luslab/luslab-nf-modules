#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for raven...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {raven} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_nanopore = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_all.fastq.gz"],
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
    raven (params.modules["raven"], ch_fastq_data)

    // Collect file names and view output
    raven.out.fasta | view
    raven.out.gfa | view

    // Verify channel counts
    assert_channel_count(raven.out.fasta, "assembly fasta", 1)
    assert_channel_count(raven.out.gfa, "assembly graph", 1)
}
