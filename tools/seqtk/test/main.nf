#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for seqtk...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.seqtk_subsample_args = ''
params.seqtk_subsample_seed = 999
params.seqtk_subsample_number = 100

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {seqtk_subsample} from '../main.nf'
include {seqtk_subsample as seqtk_subsample_pe} from '../main.nf'
include {seqtk_subseq} from '../main.nf'
include {seqtk_subseq as seqtk_subseq2} from '../main.nf'
include {decompress} from '../../../tools/luslab_linux_tools/main.nf'
include {region2bed} from '../../../tools/luslab_genome_tools/main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testDataSingleEnd = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz"]
] 

testDataPairedEnd = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile1_r2.fq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile2_r2.fq.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile3_r2.fq.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile4_r2.fq.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile5_r2.fq.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile6_r2.fq.gz"]
]

Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_fastq_single_end}

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

Channel
    .value([[:],"$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"])
    .set {ch_fasta}

Channel
    .value("$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz")
    .set {ch_fastq}

Channel
    .value("21:40000000-40100000")
    .set {ch_region} 

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run seqtk single end
    seqtk_subsample ( params.modules['seqtk_subsample'], ch_fastq_single_end )
    seqtk_subsample.out.fastq | view

    // Run seqtk paired end
    seqtk_subsample_pe ( params.modules['seqtk_subsample'], ch_fastq_paired_end )
    seqtk_subsample_pe.out.fastq | view

    // Run fasta sub sample
    region2bed ( ch_region )
    decompress( ch_fasta )
    seqtk_subseq( params.modules['seqtk_subseq'], decompress.out.fileNoMeta, region2bed.out.bed )
    seqtk_subseq.out.subset | view

    // Run fastq subsample
    seqtk_subseq2( params.modules['seqtk_subseq'], ch_fastq, region2bed.out.bed )
    seqtk_subseq2.out.subset | view

    assert_channel_count( seqtk_subsample.out.fastq, "fastq", 6)
    assert_channel_count( seqtk_subsample_pe.out.fastq, "fastq", 6)
    assert_channel_count( seqtk_subseq.out.subset, "fastq", 1)
    assert_channel_count( seqtk_subseq2.out.subset, "fastq", 1)
}