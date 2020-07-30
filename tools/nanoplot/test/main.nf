#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Log out
log.info ("Starting tests for NanoPlot...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {nanoplot} from "../main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataNanopore= [
    ["test-sample", "$baseDir/input/test-nanopore.fastq.gz"],
]

//Define test data input channels
Channel
    .from(testDataNanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run NanoPlot
    nanoplot ( ch_fastq_data )

    // Collect file names and view output
    nanoplot.out.nanoplotOutputs | view
}
