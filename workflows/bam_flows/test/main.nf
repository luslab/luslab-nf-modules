#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bam_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include { sort_index_bam } from '../../../workflows/bam_flows/main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

bam_test = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam"]
]

// Define BAM channel
Channel
    .from(bam_test)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bam}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run paired_bam_to_bedgraph
    sort_index_bam ( ch_bam )

    sort_index_bam.out.bam_bai | view

    // Check count
    assert_channel_count( sort_index_bam.out.bam_bai, "bam_bai", 2)
}