#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for raven...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {raven} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"
include {assert_line_number} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_nanopore = [
    [[sample_id:"test_sample"], "$baseDir/../../../test_data/lambda1000a/lambda_all.fastq.gz"],
]

// Define test data input channels

// Nanopore data
Channel
    .from(test_data_nanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_data}

// Line count assertion to check that automated processes produce the expected outputs
output_line_counts = [
    fasta : 2,
    gfa : 1
]

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run flye
    raven (params.modules["raven"], ch_fastq_data)

    // Collect file names and view output
    raven.out.fasta | view
    raven.out.gfa | view

    // Verify channel counts
    assert_channel_count(raven.out.fasta, "assembly fasta", 1)
    assert_channel_count(raven.out.gfa, "assembly graph", 1)

    // Verify output line count_lines.
    // First, remove the metadata from the output channel meta:output tuples,
    // renaming the meta 'key' to match output_line_counts.
    fasta_channel = raven.out.fasta.map { row -> ['fasta', row[1]] }
    gfa_channel = raven.out.gfa.map { row -> ['gfa', row[1]] }
    simplified_output_channel = fasta_channel.concat(gfa_channel)
        .view()

    assert_line_number(simplified_output_channel, "concatenated outputs", output_line_counts)
}
