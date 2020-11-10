#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for minimap2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {minimap2_index} from "../main.nf"
include {minimap2_paf} from "../main.nf"
include {minimap2_sam} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data
test_fasta = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]

test_fastq = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq.gz"],
]

// Define test data input channel
Channel
    .from(test_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_fasta}

Channel
    .from(test_fastq)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_fastq}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minimap2
    minimap2_index(params.modules["minimap2_index"], ch_test_fasta)
    minimap2_paf(params.modules["minimap2_paf"], ch_test_fasta, ch_test_fastq)
    minimap2_sam(params.modules["minimap2_sam"], ch_test_fasta, ch_test_fastq)

    // Collect and view output
    minimap2_index.out.mmi | view
    minimap2_paf.out.paf | view
    minimap2_sam.out.sam | view

    assert_channel_count(minimap2_index.out.mmi, "index", 1)
    assert_channel_count(minimap2_paf.out.paf, "overlaps", 1)
    assert_channel_count(minimap2_sam.out.sam, "overlaps", 1)
}
