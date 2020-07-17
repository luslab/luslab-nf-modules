#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for guppy...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.guppy_flowcell = "FLO-MIN106"
params.guppy_kit = "SQK-RAD002"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include guppy_basecaller from '../main.nf'
include guppy_qc from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 

testData = "$baseDir/input"

//Define test data input channel
Channel
    .fromPath(testData, checkIfExists: true)
    .set {ch_fast5}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run guppy_basecaller
    guppy_basecaller ( ch_fast5 )
    guppy_qc ( guppy_basecaller.out.summary )

    // Collect file names and view output
    guppy_basecaller.out.basecalledSeq | view
    guppy_qc.out.report | view
}