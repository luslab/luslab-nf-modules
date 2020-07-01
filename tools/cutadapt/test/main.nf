#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for cutadapt...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.cutadapt_args = '-a AGATCGGAAGAGC'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include cutadapt from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataSingleEnd= [
    ['Sample1', "$baseDir/input/sample1_r1.fq.gz"],
    ['Sample2', "$baseDir/input/sample2_r1.fq.gz"],
    ['Sample3', "$baseDir/input/sample3_r1.fq.gz"],
    ['Sample4', "$baseDir/input/sample4_r1.fq.gz"],
    ['Sample5', "$baseDir/input/sample5_r1.fq.gz"],
    ['Sample6', "$baseDir/input/sample6_r1.fq.gz"]
]

testDataPairedEnd= [
    ['Sample1', "$baseDir/input/sample1_r1.fq.gz", "$baseDir/input/sample1_r2.fq.gz" ],
    ['Sample2', "$baseDir/input/sample2_r1.fq.gz", "$baseDir/input/sample2_r2.fq.gz"],
    ['Sample3', "$baseDir/input/sample3_r1.fq.gz", "$baseDir/input/sample3_r2.fq.gz"],
    ['Sample4', "$baseDir/input/sample4_r1.fq.gz", "$baseDir/input/sample4_r2.fq.gz"],
    ['Sample5', "$baseDir/input/sample5_r1.fq.gz", "$baseDir/input/sample5_r2.fq.gz"],
    ['Sample6', "$baseDir/input/sample6_r1.fq.gz", "$baseDir/input/sample6_r2.fq.gz"]
]

//Define test data input channels

//Single end
Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_single_end}

//Paired-end
Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}
/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // Run cutadapt
    //cutadapt ( ch_fastq_single_end )
    cutadapt ( ch_fastq_paired_end )

    // Collect file names and view output
    cutadapt.out.trimmedReads | view

}