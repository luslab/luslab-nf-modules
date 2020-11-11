#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for EMBOSS...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {emboss_seqret} from "../main.nf"
params.modules["emboss_seqret"].input_format = "fasta"
params.modules["emboss_seqret"].output_format = "phylip"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data
test_alignment = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/insulin.afa"],
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
    emboss_seqret(params.modules["emboss_seqret"], ch_alignment)

    // Collect and view output
    emboss_seqret.out.out_seq | view

    assert_channel_count( emboss_seqret.out.out_seq, "converted_file", 1)
}
