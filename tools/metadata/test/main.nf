#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for metadata...")

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {fastq_metadata as meta_se; fastq_metadata as meta_pe;} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    meta_se("$baseDir/../../../test_data/metadata/single_end_test.csv")
    meta_pe("$baseDir/../../../test_data/metadata/paired_end_test.csv")

    meta_se.out.metadata | view
    meta_pe.out.metadata | view
}