#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test pipeline for multiqc")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include multiqc from '../main.nf'

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // Create test data channel from log files
    ch_testData = Channel.fromPath( "test" )

    // Run multiqc
    multiqc ( ch_testData )

    // Collect file names and view output
    multiqc.out..collect() | view
}