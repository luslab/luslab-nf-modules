#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting paraclu module testing")

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

params.paraclu_min_value = 10
params.paraclu_max_cluster_length = 200
params.paraclu_min_density_increase = 2

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include paraclu from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Defining input channels
--------------------------------------------------------------------------------------*/

// Defining test data
testData = [
    ['Sample1', "$baseDir/input/sample1.xl.bed.gz"],
    ['Sample2', "$baseDir/input/sample2.xl.bed.gz"],
    ['Sample3', "$baseDir/input/sample3.xl.bed.gz"],
    ['Sample4', "$baseDir/input/sample4.xl.bed.gz"],
    ['Sample5', "$baseDir/input/sample5.xl.bed.gz"],
    ['Sample6', "$baseDir/input/sample6.xl.bed.gz"]
]

// Define test data input channels
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_crosslinks}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run paraclu
    paraclu ( ch_crosslinks )

    // Collect file names and view output
    paraclu.out.peaks | view
}