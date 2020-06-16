#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting paraclu module testing")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include paraclu from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

testCrosslinks = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.xl.bed.gz"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.xl.bed.gz"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.xl.bed.gz"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.xl.bed.gz"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.xl.bed.gz"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.xl.bed.gz"]
]

  Channel
  .from(testCrosslinks)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .set {ch_test_crosslinks}

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // Run paraclu
    paraclu( ch_test_crosslinks )

    // Collect file names and view output
    paraclu.out.peaks | view
}