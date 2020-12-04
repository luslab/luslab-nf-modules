#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for purge_haplotigs...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {purge_haplotigs} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_bam = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/purge_haplotigs/cns_p_ctg.aligned.sd.bam"],
]
test_data_fasta = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/purge_haplotigs/cns_p_ctg.fasta"],
]

//Define test data input channels
Channel
    .from(test_data_bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bam}
Channel
    .from(test_data_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    purge_haplotigs(params.modules["purge_haplotigs"], ch_bam, ch_fasta )

    // Collect file names and view output
    purge_haplotigs.out.whatever | view

    // Verify channel counts
    assert_channel_count(purge_haplotigs.out.whatever, "whatever", 1)
}
