#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*======================
Test STAR mapping module
======================*/

// Log
log.info ("Starting tests for STAR mapping...")

// Define main params
params.modules['star_alignReads'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.verbose = true

// Define optional input
params.modules['star_alignReads'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
params.modules['star_alignReads'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"

// Module inclusions
include { star_alignReads as map_se } from '../main.nf'

// Define input channels

// Single-end test reads
testMetaDataSingleEnd = [
  [[sample_id:'Sample1'], "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep1.Unmapped.fq"],
  [[sample_id:'Sample2'], "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep2.Unmapped.fq"]
]

// Channel for single-end reads 
Channel
    .from( testMetaDataSingleEnd )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testData_single_end }

// Channel for genome index
Channel
    .value(file("$baseDir/../../../test_data/star_index/hs_chr6_1Mb/2.7.5a"))
    .set { ch_test_index_file }

// Run tests
workflow {
    // Run single-end read mapping
    log.info ("Run single-end read mapping...")
    map_se ( params.modules['star_alignReads'], ch_testData_single_end, ch_test_index_file )
}