#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl = 2

// Log out
log.info ("Starting tests for racon...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {racon} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_racon= [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq.gz"],
]

//Define test data input channels
Channel
    .from(test_data_racon)
    .map { row -> [ row[0], file(row[1], checkIfExists: true )] }
    .set {ch_reads}

Channel
    .fromPath("$baseDir/../../../test_data/lambda1000a/lambda_top10_overlaps.paf")
		.set{ch_paf_overlaps}

Channel
    .fromPath("$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta")
		.set{ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run racon on the test set
    racon(params.modules['racon'], ch_reads, ch_paf_overlaps, ch_fasta)

    // Collect file names and view output
    racon.out.fasta | view

		// Assert channel counts
		assert_channel_count( racon.out.fasta, "polished_assembly", 1 )
}
