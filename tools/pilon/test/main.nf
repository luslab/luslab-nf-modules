#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for pilon...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {pilon} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_fasta = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

test_bam = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_SRR5042715_downsampled.bam", "$baseDir/../../../test_data/lambda1000a/lambda_SRR5042715_downsampled.bam.bai"],
]

Channel
    .from(test_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

Channel
    .from(test_bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set {ch_bam}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    pilon(params.modules['pilon'], ch_fasta, ch_bam)

    // Collect file names and view output
    pilon.out.pilon | view
    pilon.out.fasta | view

    assert_channel_count( pilon.out.pilon, "pilon", 1)
    assert_channel_count( pilon.out.fasta, "polished fasta", 1)
}
