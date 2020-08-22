#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for umi_tools dedup...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.umitools_dedup_args = '--umi-separator=":"'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {umitools_dedup} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:"sample1"], "$baseDir/input/sample1.bam", "$baseDir/input/sample1.bai"],
    [[sample_id:"sample2"], "$baseDir/input/sample2.bam", "$baseDir/input/sample2.bai"]
]

//Define test data input channel
Channel
    .from(testData)
    .map{row -> [row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists:true)]]}
    .set{ch_bam_bai}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {

    // Run dedup
    umitools_dedup(ch_bam_bai)

    // Collect file names and view output
    umitools_dedup.out.dedupBam | view
    umitools_dedup.out.dedupBai | view

    assert_channel_count( umitools_dedup.out.dedupBam, "dedupBam", 2)
    assert_channel_count( umitools_dedup.out.dedupBai, "dedupBai", 2)

}