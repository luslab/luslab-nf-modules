#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for report_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include { bt2_parse } from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

report_meta_test = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/report_flows/bowtie2_stats_exp.txt"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/report_flows/bowtie2_stats_spike.txt"]
]

// Define report channel
Channel
    .from( report_meta_test )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_report_meta }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Parse bt2 report 
    bt2_parse ( ch_report_meta )

    //bt2_parse.out.key_values
}