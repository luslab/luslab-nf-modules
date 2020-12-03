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

include {augustus_run_included} from "../main.nf"
params.modules["augustus_run_included"].species = "saccharomyces_cerevisiae_S288C"
params.modules["augustus_run_included"].args = "--progress=true --softmasking=on"
params.modules["augustus_run_included"].publish_dir = "augustus_run_included"


include {augustus_run_custom} from "../main.nf"
params.modules["augustus_run_custom"].species_dir = "$baseDir/../../../test_data/augustus/s_cerevisiae_custom"
params.modules["augustus_run_custom"].species_name = "s_cerevisiae_custom"
params.modules["augustus_run_custom"].args = "--progress=true --softmasking=on"
params.modules["augustus_run_custom"].publish_dir = "augustus_run_custom"

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
    augustus_run_included(params.modules["augustus_run_included"], ch_fasta)
    augustus_run_custom(params.modules["augustus_run_custom"], ch_fasta)

    // Collect file names and view output
    augustus_run_included.out.gff | view
    augustus_run_custom.out.gff | view

    assert_channel_count(augustus_run_included.out.gff, "augustus_included_model_out", 1)
    assert_channel_count(augustus_run_custom.out.gff, "augustus_custom_model_out", 1)
}
