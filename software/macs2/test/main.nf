#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for macs2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
//params.modules['bowtie2_align'].args = '--very-sensitive'

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
// For building test dataset
//include {bowtie2_align; bowtie2_build} from '../../../tools/bowtie2/main.nf'

// Main tests
include {macs2_callpeaks} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

/*test_gen_pe_fastq = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/SRR1822153_1.fastq.gz", "$baseDir/../../../test_data/atac-seq/SRR1822153_2.fastq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/SRR1822154_1.fastq.gz", "$baseDir/../../../test_data/atac-seq/SRR1822154_2.fastq.gz"]
]

Channel
    .from(test_gen_pe_fastq)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_pe}

Channel
    .fromPath("$baseDir/../../../test_data/atac-seq/genome.fa")
    .set {ch_genome}*/

bam_input = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam", "$baseDir/../../../test_data/atac-seq/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam", "$baseDir/../../../test_data/atac-seq/sample2.bam.bai"]
]

Channel
    .from(bam_input)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_bam}

//------------------------------------------------------------------------------------

workflow {
    // Test data generation
    //bowtie2_build( params.modules['bowtie2_build'], ch_genome )
    //bowtie2_align( params.modules['bowtie2_align'], ch_fastq_pe, bowtie2_build.out.bowtieIndex.collect() )
    //bowtie2_align.out.bam | view

    // Main tests
    macs2_callpeaks( params.modules['macs2'], ch_bam )
    macs2_callpeaks.out.peaks | view
    macs2_callpeaks.out.xls | view

    assert_channel_count( macs2_callpeaks.out.peaks, "peaks", 2)
    assert_channel_count( macs2_callpeaks.out.xls, "xls", 2)
}