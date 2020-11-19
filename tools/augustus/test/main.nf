#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for AUGUSTUS...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {augustus} from "../main.nf"
params.modules["augustus"].species = "saccharomyces_cerevisiae_S288C"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_fasta = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/fasta/S_cerevisiae_chrI.fa"],
]

Channel
    .from(test_data_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    augustus(params.modules['augustus'], ch_fasta)

    // Collect file names and view output
    augustus.out.augustus | view

    assert_channel_count(augustus.out.augustus, "augustus", 1)
}
