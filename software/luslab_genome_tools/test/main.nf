#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for genome_tools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {region2bed} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Create channels of test data 
Channel
    .from("chrX:15560138-15602945")
    .set {ch_region} 

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    region2bed ( ch_region )

    assert_channel_count( region2bed.out.bed, "bed", 1)
}