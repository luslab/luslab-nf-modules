#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for get_crosslinks...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.fai = "$baseDir/input/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa.fai"
params.get_crosslinks_args = ''
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include getcrosslinks from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testMetaDataBam = [
    ['Sample1', "$baseDir/input/sample1.bam"],
    ['Sample2', "$baseDir/input/sample2.bam"]
]

testMetaDataBamBai = [
    ['Sample1', "$baseDir/input/sample1.bam", "$baseDir/input/sample1.bam.bai"],
    ['Sample2', "$baseDir/input/sample2.bam", "$baseDir/input/sample2.bam.bai"]
]

//Define test data input channels

// Fai input channel
Channel
    .fromPath(params.fai, checkIfExists: true)
    .set {ch_test_fai}

// Bam input channel
Channel
    .from(testMetaDataBam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine( ch_test_fai )
    .set {ch_test_meta_bam}

// Bam/bai input channel
Channel
    .from(testMetaDataBamBai)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists:true)] ] }
    .combine( ch_test_fai )
    .set {ch_test_meta_bam_bai}


/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run getcrosslinks
    getcrosslinks ( ch_test_meta_bam )
    //getcrosslinks ( ch_test_meta_bam_bai )

    // Collect file names and view output
    getcrosslinks.out.crosslinkBed | view
}