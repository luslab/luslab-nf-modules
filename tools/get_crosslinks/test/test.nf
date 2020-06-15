#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting getcrosslinks module testing")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include getcrosslinks from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

params.fai = [["$baseDir/input/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa.fai"]]

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testMetaDataBam = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.dedup.bam"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.dedup.bam"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.dedup.bam"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.dedup.bam"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.dedup.bam"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.dedup.bam"]
]

// Create channels of test data 

// Fai input channel
  Channel
  .from(params.fai)
  .map { row -> file(row[0], checkIfExists: true)}
  .set {ch_test_fai}

//Bam input channel
  Channel
  .from(testMetaDataBam)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .combine( ch_test_fai )
  .set {ch_test_meta_bam}

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {

    // Run getcrosslinks
    getcrosslinks( ch_test_meta_bam)

    // Collect file names and view output
    getcrosslinks.out.crosslinkBed | view
}