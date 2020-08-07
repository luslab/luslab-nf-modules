#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for htseq...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules['htseq_count'].args = '-f bam -s no -m union'

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {htseq_count} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataPairedEnd= [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/htseq/s1a_sub.fastq.gz", "$baseDir/../../../test_data/htseq/s1b_sub.fastq.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/htseq/s2a_sub.fastq.gz", "$baseDir/../../../test_data/htseq/s2b_sub.fastq.gz"]
]

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

 Channel
    .fromPath("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_gtf}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    htseq_count( params.modules['htseq_count'], ch_fastq_paired_end, ch_gtf )

    //Check count
    assert_channel_count( htseq_count.out.counts, "counts", 2)
}