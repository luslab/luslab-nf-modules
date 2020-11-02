#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for CD-HIT...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {cdhit_prot} from "../main.nf"
include {cdhit_nucl} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_fasta_nucl = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

test_fasta_prot = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

Channel
    .from(test_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    minionqc(params.modules['minionqc'], ch_sequencing_summary)

    // Collect file names and view output
    minionqc.out.minionqc_output_dir | view

    assert_channel_count( minionqc.out.minionqc_output_dir, "reads", 1)
}
