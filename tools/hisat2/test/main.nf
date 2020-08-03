#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for hisat2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {hisat2} from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/


testDataSingleEnd= [
    ['Sample1', "$baseDir/test_data/hisat2/se_1.trimmed.fq.gz"],
    ['Sample2', "$baseDir/test_data/hisat2/se_2.trimmed.fq.gz"]
]

testDataPairedEnd= [
    ['Sample1', "$baseDir/test_data/hisat2/pe_1a.trimmed.fq.gz", "$baseDir/test_data/hisat2/pe_1b.trimmed.fq.gz"],
    ['Sample2', "$baseDir/test_data/hisat2/pe_2a.trimmed.fq.gz", "$baseDir/test_data/hisat2/pe_2b.trimmed.fq.gz"],
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
    hisat2 ( ch_fastq_paired_end )

    // Collect file names and view output
    cutadapt.out.sam | view
}