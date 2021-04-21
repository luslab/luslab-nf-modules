#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for flye...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules["flye"].genome_size = "50000"
// The following flag is required for low-coverage assemblies - like the test set
// params.modules["flye"].args = "--asm-coverage 4"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {flye} from "../main.nf"
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
    flye ( params.modules["flye"], ch_fastq_data )

    // Collect file names and view output
    flye.out.fasta | view

    // Verify channel counts
    assert_channel_count(flye.out.fasta, "assembly", 1)
}
