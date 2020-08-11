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
include { star_alignReads as map_pe } from '../main.nf'

// Define input channels

// Paired-end test reads
testMetaDataPairedEnd = [
  [[sample_id:'Sample1'], "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"]
]

// Channel for genome index
Channel
    .value(file("$baseDir/../../../test_data/star_index/hs_chr6_1Mb/2.7.5a"))
    .set { ch_test_index_file }
       
// Channel for paired-end reads
Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .set { ch_testData_paired_end }

// Run tests
workflow {
    // Run paired-end read mapping
    log.info ("Run paired-end read mapping...")
    map_pe ( params.modules['star_alignReads'], ch_testData_paired_end, ch_test_index_file )
}

