#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for get_crosslinks...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {getcrosslinks as getcrosslinks_bam; getcrosslinks as getcrosslinks_bambai} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
// testDataBam = [
//     [[sample_id:"sample1"], "$baseDir/../../../test_data/bam_bai/sample1.bam"],
//     [[sample_id:"sample2"], "$baseDir/../../../test_data/bam_bai/sample2.bam"]
// ]

testDataBamBai = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bam_bai/sample1.bam", "$baseDir/../../../test_data/bam_bai/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bam_bai/sample2.bam", "$baseDir/../../../test_data/bam_bai/sample2.bam.bai"]
]

//Define test data input channels

// Fai input channel
Channel
    .value("$baseDir/../../../test_data/fai/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa.fai")
    .set {ch_fai}

// Bam input channel
// Channel
//     .from(testDataBam)
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set {ch_bam}

// Bam/bai input channel
Channel
    .from(testDataBamBai)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists:true)] ] }
    .set {ch_bam_bai}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // getcrosslinks_bam (params.modules['get_crosslinks'], ch_bam, ch_fai )
    getcrosslinks_bambai (params.modules['get_crosslinks'], ch_bam_bai, ch_fai )

    // getcrosslinks_bam.out.crosslinkBed | view
    getcrosslinks_bambai.out.crosslinkBed | view

    // assert_channel_count( getcrosslinks_bam.out.crosslinkBed, "crosslinkBed", 2)
    assert_channel_count( getcrosslinks_bambai.out.crosslinkBed, "crosslinkBed", 2)
}