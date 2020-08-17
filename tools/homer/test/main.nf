#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for homer...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.verbose = true
params.gtf = "$baseDir/../../../test_data/homer/chr1.gtf"
params.fasta = "$baseDir/../../../test_data/homer/Gallus_gallus.sub.fa"

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {homer_annotate_peaks} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
homerData = [
    [[sample_id:'S1'], "$baseDir/../../../test_data/homer/testPeaks.bed"]
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
    homer_annotate_peaks(params.modules['homer_annotate_peaks'], ch_homerData, params.fasta, params.gtf)

    // View outputs
    homer_annotate_peaks.out | view
}