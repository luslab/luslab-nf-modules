#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test pipeline for nf-core peka module")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include peka from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

testPeaks = [
  ['Sample1', "$baseDir/input/K562-TIA1-chr20.xl_peaks.bed.gz"]
]

testXls = [
  ['Sample1', "$baseDir/input/K562-TIA1-chr20.xl.bed.gz"]
]

testGenome = [["$baseDir/input/chr20.fa"]]
testGenomeIndex = [["$baseDir/input/chr20.fa.fai"]]
testRegions = [["$baseDir/input/regions_GENCODE_v30.gtf.gz"]]

Channel
  .from(testRegions)
  .map { row -> file(row[0], checkIfExists: true)}
  .set {ch_regions}

Channel
  .from(testGenomeIndex)
  .map { row -> file(row[0], checkIfExists: true)}
  .set {ch_genome_index}

Channel
  .from(testGenome)
  .map { row -> file(row[0], checkIfExists: true)}
  .set {ch_genome}

Channel
  .from(testXls)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .set {ch_bed}

Channel
  .from(testPeaks)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .join(ch_bed)
  .combine(ch_genome)
  .combine(ch_genome_index)
  .combine(ch_regions)
  .set {ch_combined}
  
/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    //Run peka
    peka( ch_combined )

    // Collect file names and view output
    peka.out.results | view
    
}