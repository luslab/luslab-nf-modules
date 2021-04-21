#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for metadata...")

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {fastq_metadata as meta_se; fastq_metadata as meta_pe; smartseq2_fastq_metadata; tenx_fastq_metadata; bam_metadata} from '../main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    meta_se("$baseDir/../../../test_data/metadata/single_end_test.csv")
    meta_pe("$baseDir/../../../test_data/metadata/paired_end_test.csv")
    smartseq2_fastq_metadata("$baseDir/../../../test_data/metadata/smartseq_test.csv")
    tenx_fastq_metadata("$baseDir/../../../test_data/metadata/10x_test.csv")
    bam_metadata("$baseDir/../../../test_data/metadata/bam_test.csv")

    meta_se.out.metadata | view
    meta_pe.out.metadata | view
    smartseq2_fastq_metadata.out.metadata | view
    tenx_fastq_metadata.out.metadata | view
    bam_metadata.out.metadata | view

    // Check count
    assert_channel_count( meta_se.out.metadata, "metadata_se", 3)
    assert_channel_count( meta_pe.out.metadata, "metadata_pe", 3)
    assert_channel_count( smartseq2_fastq_metadata.out.metadata, "metadata_ss2", 6)
    assert_channel_count( tenx_fastq_metadata.out.metadata, "metadata_10x", 2)
    assert_channel_count( bam_metadata.out.metadata, "metadata_bam", 3)
}