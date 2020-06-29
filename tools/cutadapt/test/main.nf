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

testData = [
    ['Sample1', "$baseDir/input/sample1.fq.gz"],
    ['Sample2', "$baseDir/input/sample2.fq.gz"],
    ['Sample3', "$baseDir/input/sample3.fq.gz"],
    ['Sample4', "$baseDir/input/sample4.fq.gz"],
    ['Sample5', "$baseDir/input/sample5.fq.gz"],
    ['Sample6', "$baseDir/input/sample6.fq.gz"]
]

//Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // Run cutadapt
    cutadapt ( ch_fastq )

    // Collect file names and view output
    cutadapt.out.trimmedReads | view

}