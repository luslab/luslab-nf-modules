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
    [[sample_id:"sample1"], "$baseDir/../../../test_data/umi_tools/sample1.bam", "$baseDir/../../../test_data/umi_tools/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/umi_tools/sample2.bam", "$baseDir/../../../test_data/umi_tools/sample2.bam.bai"]
]

//Define test data input channel
Channel
    .from(testData)
    .map{row -> [row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists:true)]}
    .set{ch_bam_bai}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {

    // Run dedup
    umitools_dedup(params.modules['umi_tools'], ch_bam_bai)

    // Collect file names and view output
    umitools_dedup.out.dedupBam | view

    assert_channel_count( umitools_dedup.out.dedupBam, "dedupBam", 2)

}