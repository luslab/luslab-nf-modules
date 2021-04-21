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


test_fasta_prot = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/fasta/insulin.faa"],
]

test_fasta_nucl = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/fasta/insulin.fna"],
]

Channel
    .from(test_fasta_prot)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta_prot}

Channel
    .from(test_fasta_nucl)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta_nucl}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    cdhit_prot(params.modules['cdhit_prot'], ch_fasta_prot)
    cdhit_nucl(params.modules['cdhit_nucl'], ch_fasta_nucl)

    // Collect file names and view output
    cdhit_prot.out.fasta | view
    cdhit_nucl.out.fasta | view

    assert_channel_count( cdhit_prot.out.fasta, "cdhit_prot", 1)
    assert_channel_count( cdhit_nucl.out.fasta, "cdhit_nucl", 1)
}
