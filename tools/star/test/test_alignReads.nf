#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*======================
Test STAR mapping module
======================*/

// Log
log.info ("Starting tests for STAR mapping...")

// Define main params
params.modules['star_alignReads'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM' //--quantMode TranscriptomeSAM GeneCounts' // '--quantMode GeneCounts'
params.verbose = true

// Define optional input
params.modules['star_alignReads'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
params.modules['star_alignReads'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"

// Module inclusions
include { star_alignReads as map_se } from '../main.nf'
include { star_alignReads as map_pe } from '../main.nf'

// Define input channels

// Single-end test reads
testMetaDataSingleEnd = [
  [[sample_id:'Sample1'], "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep1.Unmapped.fq"],
  [[sample_id:'Sample2'], "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep2.Unmapped.fq"]
]

// Paired-end test reads
/*testMetaDataPairedEnd = [
  [[sample_id:'Sample1'], "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"]
]*/

// Channel for single-end reads 
Channel
    .from( testMetaDataSingleEnd )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testData_single_end }

// Channel for genome index
Channel
    .value(file("$baseDir/../../../test_data/star_index/hs_chr6_1Mb/2.7.5a"))
    .set { ch_test_index_file }
       
// Channel for paired-end reads
/*Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .set { ch_testData_paired_end }*/

// Run tests
workflow {
    // Run single-end read mapping
    log.info ("Run single-end read mapping...")
    map_se ( params.modules['star_alignReads'], ch_testData_single_end, ch_test_index_file )

    // Collect file names and view output for single-end read mapping
    /*map_se.out.bamFiles.collect() | view
    map_se.out.samFiles.collect() | view
    map_se.out.sjFiles.collect() | view
    map_se.out.finalLogFiles.collect() | view
    map_se.out.outLogFiles.collect() | view
    map_se.out.progressLogFiles.collect() | view */

    // Run paired-end read mapping
    //log.info ("Run paired-end read mapping...")
    //map_pe ( params.modules['star_alignReads'], ch_testData_paired_end, ch_test_index_file )

    // Collect file names and view output for paired-end read mapping
    /*map_pe.out.bamFiles.collect() | view
    map_pe.out.samFiles.collect() | view
    map_pe.out.sjFiles.collect() | view
    map_pe.out.finalLogFiles.collect() | view
    map_pe.out.outLogFiles.collect() | view
    map_pe.out.progressLogFiles.collect() | view 
    map_pe.out.readsPerCount.collect() | view */
}

