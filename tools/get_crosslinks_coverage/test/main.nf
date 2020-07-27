#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for get_crosslinks_coverage...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include getcrosslinkscoverage from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testData = [
    ['Sample1', "$baseDir/input/sample1.xl.bed.gz"],
    ['Sample2', "$baseDir/input/sample2.xl.bed.gz"],
    ['Sample3', "$baseDir/input/sample3.xl.bed.gz"],
    ['Sample4', "$baseDir/input/sample4.xl.bed.gz"],
    ['Sample5', "$baseDir/input/sample5.xl.bed.gz"],
    ['Sample6', "$baseDir/input/sample6.xl.bed.gz"]
]

// Create channels of test data 
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bed} 

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    // Run getcrosslinkcoverage
    getcrosslinkscoverage ( ch_bed )

    // Collect file names and view output
    getcrosslinkscoverage.out.bedGraph | view
    getcrosslinkscoverage.out.normBedGraph | view
}