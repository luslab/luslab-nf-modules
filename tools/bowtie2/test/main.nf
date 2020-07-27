#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bowtie2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.bowtie2_args = '--very-sensitive'
params.bowtie2_build_args = ''

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {bowtie2_align; bowtie2_build} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataPairedEnd= [
    ['sample1', "$baseDir/../../../test_data/fastq/ENCFF038BYR.sub.fastq.gz", "$baseDir/../../../test_data/fastq/ENCFF721JZG.sub.fastq.gz"]
]

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

 Channel
    .from("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    bowtie2_build( ch_fasta )

    bowtie2_align ( ch_fastq_paired_end.combine(bowtie2_build.out.bowtieIndex.map{x -> [x]}) )

    bowtie2_align.out.alignedReads | view
}

    //ch_fastq_paired_end.merge(bowtie2_build.out.bowtieIndex.map{x -> [x]})
      //  .subscribe {log.info("$it")}