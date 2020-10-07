#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for SEACR...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules['seacr'].args = "norm stringent"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include { seacr } from '../main.nf'

 
/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
hist_test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bedgraph"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bedgraph"]
]

// control_test_data = [
//     ["$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bedgraph"],
//     ["$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bedgraph"]
// ]

control_test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bedgraph"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bedgraph"]
]

// Define Channel
Channel
    .from( hist_test_data )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_hist }

Channel
    .from( control_test_data )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_control }

/*------------------------------------------------------------------------------------*/
/* Run tests
/*------------------------------------------------------------------------------------*/

workflow { 
    // Run SEACR module
    seacr ( params.modules['seacr'], ch_hist, ch_control)
}
