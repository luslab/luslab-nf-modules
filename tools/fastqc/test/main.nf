#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for fastqc...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.fastqc_args = ''
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include fastqc from '../main.nf' addParams(fastqc_reportname: 'test')

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
    // Run fastqc
    //fastqc ( ch_fastq_single_end )
    fastqc ( ch_fastq_paired_end )


    // Collect file names and view output
    fastqc.out.report | view
}