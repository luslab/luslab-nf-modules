#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting icount trimming test pipeline")

/* Define global params
--------------------------------------------------------------------------------------*/


/* Module inclusions 
--------------------------------------------------------------------------------------*/

include icount from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testMetaData = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.xl.bed.gz"]//,
//  ['Sample 2', "$baseDir/input/prpf8_ctrl_rep2.xl.bed.gz"],
//  ['Sample 3', "$baseDir/input/prpf8_ctrl_rep4.xl.bed.gz"],
//  ['Sample 4', "$baseDir/input/prpf8_eif4a3_rep1.xl.bed.gz"],
//  ['Sample 5', "$baseDir/input/prpf8_eif4a3_rep2.xl.bed.gz"],
//  ['Sample 6', "$baseDir/input/prpf8_eif4a3_rep4.xl.bed.gz"]
]

testSegPath = [
  ["$baseDir/input/segmentation.gtf.gz"]
]

// Create channels of test data 
Channel
  .from(testSegPath)
  .map { row -> file(row[0], checkIfExists: true) }
  .set {ch_test_seg}

  Channel
  .from(testMetaData)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .combine(ch_test_seg)
  .set {ch_test_meta}

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    icount( ch_test_meta ) 

    icount.out.peaks | view
    icount.out.peak_scores | view
    icount.out.clusters | view
}