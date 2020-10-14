#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bowtie2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules['bowtie2_align'].args = '--very-sensitive'
params.modules['bowtie2_align'].summary_name = "bowtie_summary"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {bowtie2_align; bowtie2_build} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataPairedEnd= [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fastq/ENCFF038BYR.sub.fastq.gz", "$baseDir/../../../test_data/fastq/ENCFF721JZG.sub.fastq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/fastq/ENCFF038BYR.sub.fastq.gz", "$baseDir/../../../test_data/fastq/ENCFF721JZG.sub.fastq.gz"]
]

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

 Channel
    .fromPath("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    bowtie2_build( params.modules['bowtie2_build'], ch_fasta )

    bowtie2_align( params.modules['bowtie2_align'], ch_fastq_paired_end, bowtie2_build.out.bowtieIndex.collect() )

    bowtie2_align.out.bam | view
    bowtie2_align.out.report | view

    //Check count
    assert_channel_count( bowtie2_align.out.sam, "sam", 0)
    assert_channel_count( bowtie2_align.out.bam, "bam", 2)
    assert_channel_count( bowtie2_align.out.unmapped_fq_pe, "unmapped_fq_pe", 0)
    assert_channel_count( bowtie2_align.out.unmapped_fq_s, "unmapped_fq_s", 0)
}

    //ch_fastq_paired_end.merge(bowtie2_build.out.bowtieIndex.map{x -> [x]})
      //  .subscribe {log.info("$it")}