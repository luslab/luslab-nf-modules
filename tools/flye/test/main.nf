#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Log out
log.info ("Starting tests for Flye...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules['flye'].genome_size = 200000

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {flye} from "../main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataNanopore= [
    [[sample_id:"test-sample"], "/Users/alex/dev/repos/luslab-nf-modules/test_data/flye/fastq_runid_b002751b18e298acfe7b2ec51dfaa0961b5d290e_0_0.sub.fastq.gz"],
]

//Define test data input channels

//Single end
Channel
    .from(testDataNanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_data}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    flye ( params.modules['flye'], ch_fastq_data )

    // Collect file names and view output
    flye.out.assemblyFasta | view
}
