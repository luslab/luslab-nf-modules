#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for BLAST...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {blast_makeblastdb} from "../main.nf"
params.modules["blast_makeblastdb"].dbtype = "nucl"
include {blast_blastn} from "../main.nf"
include {blast_windowmasker} from "../main.nf"
include {blast_asn_to_tab} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
test_reference_genome = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]
test_query_sequence = [
    [[sample_id:"query_sequence"], "$baseDir/../../../test_data/lambda1000a/lambda_subsequence.fasta"],
]

// Define test data input channels
Channel
    .from(test_reference_genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_ref_genome}

Channel
    .from(test_query_sequence)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_sequence}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Work
    blast_makeblastdb(params.modules['blast_makeblastdb'], ch_ref_genome)
    blast_blastn(params.modules['blast_blastn'], ch_ref_genome, blast_makeblastdb.out.blast_db, ch_query_sequence)
    blast_asn_to_tab(params.modules['blast_asn_to_tab'], blast_blastn.out.blast_output, blast_makeblastdb.out.blast_db)
    blast_windowmasker(params.modules['blast_windowmasker'], ch_ref_genome)

    // Print outputs
    blast_makeblastdb.out.blast_db | view
    blast_blastn.out.blast_output | view
    blast_asn_to_tab.out.blast_output | view
    blast_windowmasker.out.fasta | view

    // Assert channel counts
    assert_channel_count(blast_makeblastdb.out.blast_db, "blast_db", 1)
    assert_channel_count(blast_blastn.out.blast_output, "initial_blast", 1)
    assert_channel_count(blast_asn_to_tab.out.blast_output, "reformatted_blast", 1)
    assert_channel_count(blast_windowmasker.out.wmstat, "windowmasker", 1)

}
