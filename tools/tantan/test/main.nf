#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for tantan...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {tantan; tantan_to_GFF3} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testData = [
    [ [sample_id:"sample1"]
    , "$baseDir/../../../test_data/homer/Gallus_gallus.sub.fa"]
]

// Define test data input channel
channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    tantan (params.modules['tantan'], ch_test_fasta)
    tantan_to_GFF3(params.modules['tantan_to_GFF3'], tantan.out.tantanRepeats)

    tantan.out.tantanRepeats | view
    tantan_to_GFF3.out.tantanRepeats | view

    //Check counts
    assert_channel_count( tantan.out.tantanRepeats, "repeats", 1)
    assert_channel_count( tantan_to_GFF3.out.tantanRepeats, "gff3 repeats", 1)
}
