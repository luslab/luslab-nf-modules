#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for last...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {last_make_index} from "../main.nf"
include {last_train_reads_on_genome} from "../main.nf"
include {last_align_reads_to_genome} from "../main.nf"
include {last_convert_maf_to_sam} from "../main.nf"
include {last_make_dotplot} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/last/E_coli_K-12.fna"],
]

// Define regions file input channel
Channel
    .fromPath("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_test_fasta}

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_fastq}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run last
    bedtools_intersect (params.modules['last'], ch_test_fastq, ch_test_fasta )

    // Collect file names and view output
    last.out.mappedReads | view

    // Double check the channel count
    assert_channel_count( busco_genome.out.report, "busco", 1 )
}
