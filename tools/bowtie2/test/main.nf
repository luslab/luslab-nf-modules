#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for bowtie2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.bowtie_args = ''
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {bowtie2_align; bowtie2_build} from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testDataPairedEnd= [
    ['Sample1', "$baseDir/input/sample1_r1.fq.gz", "$baseDir/input/sample1_r2.fq.gz"]
]

//Channel
//    .from(testDataPairedEnd)
//    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
//    .set {ch_fastq_paired_end}

Channel
    .from("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    bowtie2_build( ch_fasta )

    // Run cutadapt
    //cutadapt ( ch_fastq_single_end )
   // bowtie2_align ( ch_fastq_paired_end )

    // Collect file names and view output
   // cutadapt.out.trimmedReads | view
}