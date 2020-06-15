#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting Bedtools intersect module test")

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include bedtools_intersect from '../main.nf'
 
/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

params.regions_file = [["$baseDir/input/regions_GENCODE_v30.gtf.gz"]]

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

testMetaCrosslinks = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.xl.bed.gz"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.xl.bed.gz"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.xl.bed.gz"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.xl.bed.gz"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.xl.bed.gz"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.xl.bed.gz"]
]

  // Regions file channel
  Channel
  .from(params.regions_file)
  .map { row -> file(row[0], checkIfExists: true)}
  .set {ch_test_regions_file}

  //Bedtools input channel
  Channel
  .from(testMetaCrosslinks)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .combine(ch_test_regions_file)
  .set {ch_test_meta_crosslinks}

/*------------------------------------------------------------------------------------*/
// Run workflow 

workflow {
    // Run bedtools_intersect
    bedtools_intersect( ch_test_meta_crosslinks )

    // Collect file names and view output
    bedtools_intersect.out.annotatedBed | view
}

