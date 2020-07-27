#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for iCount...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.seg = "$baseDir/input/segmentation.gtf.gz"
params.icount_half_window = '3'
params.icount_fdr = 0.05

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include icount from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
  ['Sample1', "$baseDir/input/sample1.xl.bed.gz"],
  ['Sample 2', "$baseDir/input/sample2.xl.bed.gz"]//,
//  ['Sample 3', "$baseDir/input/prpf8_ctrl_rep4.xl.bed.gz"],
//  ['Sample 4', "$baseDir/input/prpf8_eif4a3_rep1.xl.bed.gz"],
//  ['Sample 5', "$baseDir/input/prpf8_eif4a3_rep2.xl.bed.gz"],
//  ['Sample 6', "$baseDir/input/prpf8_eif4a3_rep4.xl.bed.gz"]
]

// Define test data input channels 

// Seg file channel
Channel
    .fromPath(params.seg, , checkIfExists: true)
    .set {ch_seg}

// Bed/seg channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine(ch_seg)
    .set {ch_bed}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run iCount
    icount( ch_bed ) 

    // Collect file names and view output
    icount.out.peaks | view
    icount.out.peak_scores | view
    icount.out.clusters | view
}