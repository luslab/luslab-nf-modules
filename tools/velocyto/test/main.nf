#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for velocyto...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {velocyto_run_smartseq2} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

bam = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/cutadapt/sample1_r1.fq.gz"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/cutadapt/sample2_r1.fq.gz"],
    [[sample_id:'Sample3'], "$baseDir/../../../test_data/cutadapt/sample3_r1.fq.gz"],
    [[sample_id:'Sample4'], "$baseDir/../../../test_data/cutadapt/sample4_r1.fq.gz"],
    [[sample_id:'Sample5'], "$baseDir/../../../test_data/cutadapt/sample5_r1.fq.gz"],
    [[sample_id:'Sample6'], "$baseDir/../../../test_data/cutadapt/sample6_r1.fq.gz"]
]

//Define test data input channels

//Single end
Channel
    .from(bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bam}

Channel
    .from()

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // Run cutadapt
    //cutadapt ( ch_fastq_single_end )
    velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_bam )

    // Collect file names and view output
    velocyto_run_smartseq2.out.velocyto | view

    //Check count
    assert_channel_count( velocyto_run_smartseq2.out.velocyto, "x", 6)
}