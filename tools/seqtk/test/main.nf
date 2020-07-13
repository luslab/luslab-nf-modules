#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

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

include seqtk_subsample from '../main.nf'
include seqtk_subsample as seqtk_subsample_pe from '../main.nf'
include seqtk_subseq from '../main.nf'
include decompress_noid from '../../../tools/luslab_file_tools/main.nf'
include region2bed from '../../../tools/luslab_genome_tools/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testDataSingleEnd = [
    ['Sample1', "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz"],
    ['Sample2', "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz"],
    ['Sample3', "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz"],
    ['Sample4', "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz"],
    ['Sample5', "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz"],
    ['Sample6', "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz"]
] 

testDataPairedEnd = [
    ['Sample1', "$baseDir/../../../test_data/fastq/readfile1_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile1_r2.fq.gz"],
    ['Sample2', "$baseDir/../../../test_data/fastq/readfile2_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile2_r2.fq.gz"],
    ['Sample3', "$baseDir/../../../test_data/fastq/readfile3_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile3_r2.fq.gz"],
    ['Sample4', "$baseDir/../../../test_data/fastq/readfile4_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile4_r2.fq.gz"],
    ['Sample5', "$baseDir/../../../test_data/fastq/readfile5_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile5_r2.fq.gz"],
    ['Sample6', "$baseDir/../../../test_data/fastq/readfile6_r1.fq.gz", "$baseDir/../../../test_data/fastq/readfile6_r2.fq.gz"]
]

testSubseq = []

Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_fastq_single_end}

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

Channel
    .from("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_fasta}

Channel
    .from("21:40000000-40100000")
    .set {ch_region} 

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run seqtk single end
    seqtk_subsample ( ch_fastq_single_end )
    seqtk_subsample.out.sampledReads | view

    // Run seqtk paired end
    seqtk_subsample_pe ( ch_fastq_paired_end )
    seqtk_subsample_pe.out.sampledReads | view

    // Run fasta sub sample
    region2bed ( ch_region )
    decompress_noid( ch_fasta )
    seqtk_subseq( decompress_noid.out.file.combine(region2bed.out.bedFile) )
    seqtk_subseq.out.subsetFasta | view
}