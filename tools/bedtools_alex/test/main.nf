#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bedtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

// -wa = entry in A for each overlap
// -wb = entry in B for each overlap
params.bedtools_intersect_args = '-wa -wb'
params.bedtools_subtract_args = '-A'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {bedtools_intersect; bedtools_subtract} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
chipData = [
    ['K4me1', "./input/ss8-K4me1_R1_peaks.broadPeak"],
    ['K4me3', "./input/ss8-K4me3_R1_peaks.broadPeak"],
    ['K27Ac', "./input/ss8-K27Ac_R1_peaks.broadPeak"],
    ['K27me3', "./input/ss8-K27me3_R1_peaks.broadPeak"]
]

atacData = [
    ['ATAC_R1', "./input/ss8_R1.mLb.clN_peaks.broadPeak"]
]

// Define test data input channel
Channel
    .from(chipData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_chipData}

ch_chipData
    .filter { it[0] == 'K27Ac' }
    .set {ch_K27Ac}

ch_chipData
    .filter { it[0] == 'K27me3' }
    .set {ch_K27me3}

Channel
    .from(atacData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_atacData}



/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // intersect peak files
    bedtools_intersect(ch_K27Ac, ch_atacData)

    bedtools_subtract(bedtools_intersect.out, ch_K27me3)

    // View outputs
    bedtools_subtract.out | view
}