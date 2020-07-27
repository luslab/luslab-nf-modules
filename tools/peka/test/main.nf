#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for peka...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

// Params passed as input
params.peka_testGenome = "$baseDir/input/chr20.fa"
params.peka_testGenomeIndex = "$baseDir/input/chr20.fa.fai"
params.peka_testRegions = "$baseDir/input/regions_GENCODE_v30_chr20.gtf.gz"

// Other params
params.peka_window = 40
params.peka_window_distal = 150
params.peka_kmer_length = 4
params.peka_top_n = 20
params.peka_percentile = 0.7
params.peka_clusters = 5
params.peka_smoothing = 6
params.peka_all_outputs = "False"
params.peka_regions_selection = "None"
params.peka_subsample = "True"
params.peka_repeats = "masked"

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include peka from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testPeaks = [
  ['Sample1', "$baseDir/input/K562-TIA1-chr20.xl_peaks.bed.gz"]
]

testXls = [
  ['Sample1', "$baseDir/input/K562-TIA1-chr20.xl.bed.gz"]
]

// Create channels of test data 

// Params channels
Channel
    .fromPath(params.peka_testRegions, checkIfExists: true)
    .set{ch_regions}

Channel
    .fromPath(params.peka_testGenomeIndex, checkIfExists: true)
    .set{ch_genome_index}

Channel
    .fromPath(params.peka_testGenome, checkIfExists: true)
    .set{ch_genome}

// Test data channels
Channel
    .from(testXls)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set{ch_bed}

Channel
    .from(testPeaks)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .join(ch_bed)
    .combine(ch_genome)
    .combine(ch_genome_index)
    .combine(ch_regions)
    .set{ch_combined}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    //Run peka
    peka ( ch_combined )

    // Collect file names and view output
    peka.out.results | view
}