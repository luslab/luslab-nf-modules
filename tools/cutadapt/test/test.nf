#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting Cutadapt trimming test pipeline")

/* Define global params
--------------------------------------------------------------------------------------*/

params.cutadapt_output_prefix = 'trimmed_'

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include cutadapt from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
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
    // Run cutadapt
    cutadapt(ch_test_meta)

    // Collect file names and view output
    cutadapt.out.trimmedReads | view

}