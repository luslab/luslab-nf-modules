#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for homer...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.verbose = true
params.gtf = "$baseDir/input/chr1.gtf"
params.fasta = "$baseDir/input/Gallus_gallus.sub.fa"

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {homer_annotatePeaks} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
homerData = [
    ['S1', "./input/testPeaks.bed"]
]

// Define test data input channel
Channel
    .from(homerData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_homerData}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // annotate peak files
    homer_annotatePeaks(ch_homerData, params.fasta, params.gtf)

    // View outputs
    homer_annotatePeaks.out | view
}