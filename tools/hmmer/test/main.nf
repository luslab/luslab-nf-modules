#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for HMMER...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {hmmer_hmmscan} from "../main.nf"
params.modules["hmmer_hmmscan"].db = "$baseDir/../../../test_data/hmmer/Pfam_subset.hmm"
include {hmmer_hmmsearch} from "../main.nf"
params.modules["hmmer_hmmsearch"].db = "$baseDir/../../../test_data/hmmer/Pfam_subset.hmm"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_genome = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/S_cerevisiae_prot.fa"],
]

Channel
    .from(test_data_genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run BUSCO on the test genome FASTA file
    hmmer_hmmscan( params.modules["hmmer_hmmscan"], ch_fasta )
    hmmer_hmmsearch( params.modules["hmmer_hmmsearch"], ch_fasta )

    // Confirm the outputs of the above command
    hmmer_hmmscan.out.tbl | view
    hmmer_hmmsearch.out.tbl | view

    // Double check the channel count
    assert_channel_count( hmmer_hmmscan.out.tbl, "hmmscan", 1 )
    assert_channel_count( hmmer_hmmsearch.out.tbl, "hmmscan", 1 )
}
