#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for mafft...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {mafft} from "../main.nf"
// Make sure that the alignment options passed to mafft are appropriate.
// You can add to "args" as necessary.
params.modules["mafft"].args = "--auto"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data
test_fasta = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/insulin.faa"],
]

// Define test data input channel
Channel
    .from(test_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minimap2
    mafft(params.modules["mafft"], ch_fasta )

    // Collect and view output
    mafft.out.mfa | view

    assert_channel_count( mafft.out.mfa, "alignment", 1)
}
