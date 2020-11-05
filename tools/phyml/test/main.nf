#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for phyml...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {phyml} from "../main.nf"
params.modules["phyml"].data_type = "aa"
params.modules["phyml"].args = "-b 50 -o n"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data
test_alignment = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/insulin.phylip"],
]

// Define test data input channel
Channel
    .from(test_alignment)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_alignment}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minimap2
    phyml(params.modules["phyml"], ch_alignment)

    // Collect and view output
    phyml.out.phyml | view

    assert_channel_count( phyml.out.phyml, "tree", 1)
}
