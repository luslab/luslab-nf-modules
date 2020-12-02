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

include {augustus_run} from "../main.nf"
params.modules["augustus_run"].species = "s_cerevisiae_custom"
params.modules["augustus_run"].args = "--progress=true --softmasking=on"

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
    augustus_run(params.modules["augustus_run"], ch_fasta)

    // Collect file names and view output
    augustus_run.out.gff | view

    assert_channel_count(augustus_run.out.gff, "augustus_out", 1)
}
