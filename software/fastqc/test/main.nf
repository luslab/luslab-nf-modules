#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for fastqc...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {fastqc as fastqc_single; fastqc as fastqc_paired;} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
test_data_single_end = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz"]
] 

test_data_paired_end = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile1_r2.fq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile2_r2.fq.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile3_r2.fq.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile4_r2.fq.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile5_r2.fq.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile6_r2.fq.gz"]
] 

//Define test data input channel
Channel
    .from(test_data_single_end)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_fastq_single_end}

Channel
    .from(test_data_paired_end)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    fastqc_single ( params.modules['fastqc_se'], ch_fastq_single_end )
    fastqc_paired ( params.modules['fastqc_pe'], ch_fastq_paired_end )

    fastqc_single.out.report | view
    fastqc_paired.out.report | view

    assert_channel_count( fastqc_single.out.report, "fastqc_single", 6)
    assert_channel_count( fastqc_paired.out.report, "fastqc_paired", 6)
}