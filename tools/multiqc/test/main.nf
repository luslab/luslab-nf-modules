#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting test pipeline for multiqc")

params.verbose = true

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {multiqc} from '../main.nf'

/*------------------------------------------------------------------------------------*/

ch_testData = Channel.fromPath( "$baseDir/input/*.zip" )

// Run workflow
workflow {
    // Run multiqc
    multiqc ( ch_testData.collect() )

    // Collect file names and view output
    multiqc.out.report | view
}