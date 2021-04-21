#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for iCount...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {icount} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
  [[sample_id:"sample1"], "$baseDir/../../../test_data/icount/sample1.xl.bed.gz"],
  [[sample_id:"sample2"], "$baseDir/../../../test_data/icount/sample2.xl.bed.gz"]
]

// Define test data input channels 

// Seg file channel
Channel
    .value("$baseDir/../../../test_data/icount/segmentation.gtf.gz")
    .set {ch_seg}

// Bed/seg channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bed}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run iCount
    icount( params.modules['icount'], ch_bed, ch_seg) 

    // Collect file names and view output
    icount.out.peaks | view
    icount.out.peak_scores | view
    icount.out.clusters | view

    assert_channel_count( icount.out.peaks, "peaks", 2)
    assert_channel_count( icount.out.peak_scores, "peak_scores", 2)
    assert_channel_count( icount.out.clusters, "clusters", 2)
}