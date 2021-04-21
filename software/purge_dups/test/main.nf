#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for purge_dups")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {purge_dups} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testFasta = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

testPaf = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/lambda1000a/lambda_top10_overlaps.paf"]
]

channel
    .from( testFasta )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testFasta }

channel
    .from( testPaf )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testPaf }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    purge_dups( params.modules['purge_dups'], ch_testFasta, ch_testPaf )
    assert_channel_count( purge_dups.out.purged_fasta, "split", 1 )
    assert_channel_count( purge_dups.out.hap_fasta, "split", 1 )
}
