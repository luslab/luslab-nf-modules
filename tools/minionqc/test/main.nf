#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for minionQC...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {minionqc} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_sequencing_summary = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lamda1000a/lambda_top10.sequence_summary.txt"],
]

Channel
    .from(test_data_sequencing_summary)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_sequencing_summary}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    minionqc(params.modules['minionqc'], ch_sequencing_summary)

    // Collect file names and view output
    minionqc.out.minionqc_output_dir | view

    assert_channel_count( minionqc.out.minionqc_output_dir, "reads", 1)
}
