#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting paraclu module testing")

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {paraclu} from '../main.nf'
include {assert_line_number} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Defining input channels
--------------------------------------------------------------------------------------*/

// Defining test data
testData = [
    ['Sample1', "$baseDir/../../../test_data/crosslinks/sample1.xl.bed.gz"],
    ['Sample2', "$baseDir/../../../test_data/crosslinks/sample2.xl.bed.gz"],
    ['Sample3', "$baseDir/../../../test_data/crosslinks/sample3.xl.bed.gz"],
    ['Sample4', "$baseDir/../../../test_data/crosslinks/sample4.xl.bed.gz"],
    ['Sample5', "$baseDir/../../../test_data/crosslinks/sample5.xl.bed.gz"],
    ['Sample6', "$baseDir/../../../test_data/crosslinks/sample6.xl.bed.gz"]
]

// Define test data input channels
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_crosslinks}

output_line_counts = [
    Sample1: 7,
    Sample2: 0,
    Sample3: 0,
    Sample4: 4,
    Sample5: 0,
    Sample6: 0
]

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run paraclu
    paraclu (params.modules["paraclu"], ch_crosslinks )

    // Collect file names and view output
    paraclu.out.peaks | view

    assert_line_number( paraclu.out.peaks, "paraclu_out_peaks", output_line_counts )
}