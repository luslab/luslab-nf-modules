#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

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

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataGuppyOut= [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/guppy_sequencing_summary/sequencing_summary.txt"],
]

//Define test data input channels
Channel
    .from(testDataGuppyOut)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_guppyout_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    minionqc (params.modules['minionqc'], ch_guppyout_data )

    // Collect file names and view output
    minionqc.out.minionqcOutputs | view
}
