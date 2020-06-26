#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test pipeline for fastqc")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include fastqc from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Input files
--------------------------------------------------------------------------------------*/

testMetaData = [
  ['Sample1', "$baseDir/input/readfile1.fq.gz"],
  ['Sample2', "$baseDir/input/readfile2.fq.gz"],
  ['Sample3', "$baseDir/input/readfile3.fq.gz"],
  ['Sample4', "$baseDir/input/readfile4.fq.gz"],
  ['Sample5', "$baseDir/input/readfile5.fq.gz"],
  ['Sample6', "$baseDir/input/readfile6.fq.gz"]
] 


// Create channels of test data 
Channel
  .from(testMetaData)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .set {ch_test_meta} 

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {

    // Run fastqc
    fastqc( ch_test_meta )

    // Collect file names and view output
    fastqc.out.report | view
}