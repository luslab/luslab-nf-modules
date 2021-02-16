#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2
// Log
log.info ("Starting tests for test_flows...")

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {assert_channel_count as pass_test; assert_channel_count as fail_test} from '../main.nf'
include {assert_line_number} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .from(1,2,3,4,5,6,7,8,9,10)
    .set {ch_items1}

line_count_test_data = [
    ['yeast', "$baseDir/../../../test_data/fasta/S_cerevisiae_prot.fa"],
    ['human', "$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"]
]

Channel
    .from(line_count_test_data)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_wc_test}

expected_line_counts = [
    yeast: 12004,
    human: 802166
]

//------------------------------------------------------------------------------------

workflow {
    // Run the test with a pass condition
    pass_test( ch_items1, "test_channel", 10 )

    assert_line_number(ch_wc_test, "wc_test_channel", expected_line_counts)

    // Run the test with the fail condition (DOES NOT WORK RIGHT NOW)
    // isException = false
    // try {
    //     fail_test( ch_items1, "test_channel", 9 )
    // } catch(Exception e) {
    //     isException = true
    //     log.info("here")
    // }

    //if(!isException) {
   //     exit 1, "Exception was not thrown for fail condition"
    //}
}