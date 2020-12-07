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

include {purge_haplotigs_hist} from "../main.nf"
include {purge_haplotigs_minima_2} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_bam = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/cns_p_ctg.alignedChained.0.015.bam"],
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
    purge_haplotigs_hist(params.modules["purge_haplotigs"], ch_bam, ch_fasta)
    purge_haplotigs_minima_2(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs)

    // Collect file names and view output
    //purge_haplotigs.out.purge_haplotigs | view

    // Verify channel counts
    //assert_channel_count(purge_haplotigs.out.purge_haplotigs, "whatever", 1)
}
