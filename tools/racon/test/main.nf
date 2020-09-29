#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Log out
log.info ("Starting tests for racon...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {racon} from "../main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataRacon= [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq.gz",
        "$baseDir/../../../test_data/lambda1000a/lambda_top10_overlaps.paf",
        "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

//Define test data input channels
Channel
    .from(testDataRacon)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true), file(row[3], checkIfExists: true) ] }
    .set {ch_racon_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run racon on the test set
    racon(params.modules['racon'], ch_racon_data)

    // Collect file names and view output
    racon.out.fasta | view

}
