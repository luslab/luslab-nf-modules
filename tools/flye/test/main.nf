#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Log out
log.info ("Starting tests for Flye...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {flye} from "../main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataNanopore= [
    ["test-sample", "60m", "$baseDir/input/test-nanopore.fastq.gz"],
]

//Define test data input channels

//Single end
Channel
    .from(testDataNanopore)
    .map { row -> [ row[0], row[1], file(row[2], checkIfExists: true) ] }
    .set {ch_fastq_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    flye ( ch_fastq_data )

    // Collect file names and view output
    flye.out.flyeAssembly | view
}
