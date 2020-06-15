#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting getcrosslinkscoverage module testing")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include getcrosslinkscoverage from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testMetaData = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.xl.bed.gz"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.xl.bed.gz"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.xl.bed.gz"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.xl.bed.gz"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.xl.bed.gz"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.xl.bed.gz"]
]

// Create channels of test data 
 Channel
  .from(testMetaData)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .set {ch_test_meta} 

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    // Run getcrosslinkcoverage
    getcrosslinkscoverage( ch_test_meta )

    // Collect file names and view output
    getcrosslinkscoverage.out.bedGraph | view
    getcrosslinkscoverage.out.normBedGraph | view
}