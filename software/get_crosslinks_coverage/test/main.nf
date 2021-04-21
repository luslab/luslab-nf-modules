#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for get_crosslinks_coverage...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {get_crosslinks_coverage} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

"$baseDir/../../../test_data/bed/sample1.xl.bed.gz"

//Define test data 
testData = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed/sample1.xl.bed.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed/sample2.xl.bed.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/bed/sample3.xl.bed.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/bed/sample4.xl.bed.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/bed/sample5.xl.bed.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/bed/sample6.xl.bed.gz"]
]

// Create channels of test data 
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bed} 

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    // Run getcrosslinkcoverage
    get_crosslinks_coverage ( params.modules['get_crosslinks_coverage'], ch_bed )

    // Collect file names and view output
    get_crosslinks_coverage.out.bedGraph | view
    get_crosslinks_coverage.out.normBedGraph | view

    assert_channel_count(get_crosslinks_coverage.out.bedGraph, "bedgraph", 6)
    assert_channel_count(get_crosslinks_coverage.out.normBedGraph, "norm_bedgraph", 6)
}