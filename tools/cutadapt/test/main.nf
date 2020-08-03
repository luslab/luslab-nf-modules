#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for cutadapt...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.modules['cutadapt'].args = '-a AGATCGGAAGAGC'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {cutadapt} from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataSingleEnd= [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/cutadapt/sample1_r1.fq.gz"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/cutadapt/sample2_r1.fq.gz"],
    [[sample_id:'Sample3'], "$baseDir/../../../test_data/cutadapt/sample3_r1.fq.gz"],
    [[sample_id:'Sample4'], "$baseDir/../../../test_data/cutadapt/sample4_r1.fq.gz"],
    [[sample_id:'Sample5'], "$baseDir/../../../test_data/cutadapt/sample5_r1.fq.gz"],
    [[sample_id:'Sample6'], "$baseDir/../../../test_data/cutadapt/sample6_r1.fq.gz"]
]

testDataPairedEnd= [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/cutadapt/sample1_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample1_r2.fq.gz" ],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/cutadapt/sample2_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample2_r2.fq.gz"],
    [[sample_id:'Sample3'], "$baseDir/../../../test_data/cutadapt/sample3_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample3_r2.fq.gz"],
    [[sample_id:'Sample4'], "$baseDir/../../../test_data/cutadapt/sample4_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample4_r2.fq.gz"],
    [[sample_id:'Sample5'], "$baseDir/../../../test_data/cutadapt/sample5_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample5_r2.fq.gz"],
    [[sample_id:'Sample6'], "$baseDir/../../../test_data/cutadapt/sample6_r1.fq.gz", "$baseDir/../../../test_data/cutadapt/sample6_r2.fq.gz"]
]

//Define test data input channels

//Single end
Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_single_end}

//Paired-end
Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}
/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // Run cutadapt
    //cutadapt ( ch_fastq_single_end )
    cutadapt ( params.modules['cutadapt'], ch_fastq_paired_end )

    // Collect file names and view output
    cutadapt.out.fastq | view
}