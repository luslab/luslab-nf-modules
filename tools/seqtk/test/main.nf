#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for seqtk...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.seqtk_subsample_args = ''
params.seqtk_subsample_seed = 999
params.seqtk_subsample_number = 100
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include seqtk_subsample from '../main.nf'
include seqtk_subsample as seqtk_subsample_pe from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testDataSingleEnd = [
    ['Sample1', "$baseDir/input/readfile1_r1.fq.gz"],
    ['Sample2', "$baseDir/input/readfile2_r1.fq.gz"],
    ['Sample3', "$baseDir/input/readfile3_r1.fq.gz"],
    ['Sample4', "$baseDir/input/readfile4_r1.fq.gz"],
    ['Sample5', "$baseDir/input/readfile5_r1.fq.gz"],
    ['Sample6', "$baseDir/input/readfile6_r1.fq.gz"]
] 

testDataPairedEnd = [
    ['Sample1', "$baseDir/input/readfile1_r1.fq.gz", "$baseDir/input/readfile1_r2.fq.gz"],
    ['Sample2', "$baseDir/input/readfile2_r1.fq.gz", "$baseDir/input/readfile2_r2.fq.gz"],
    ['Sample3', "$baseDir/input/readfile3_r1.fq.gz", "$baseDir/input/readfile3_r2.fq.gz"],
    ['Sample4', "$baseDir/input/readfile4_r1.fq.gz", "$baseDir/input/readfile4_r2.fq.gz"],
    ['Sample5', "$baseDir/input/readfile5_r1.fq.gz", "$baseDir/input/readfile5_r2.fq.gz"],
    ['Sample6', "$baseDir/input/readfile6_r1.fq.gz", "$baseDir/input/readfile6_r2.fq.gz"]
] 

//Define test data input channel
Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_fastq_single_end}

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

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
}