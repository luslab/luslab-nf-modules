#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*==============================
Test STAR genome indexing module
==============================*/

// Log
log.info ("Starting tests for STAR genome indexing...")

// Define main params
params.modules['star_genome_generate'].args = '--genomeSAindexNbases 9' 
params.verbose = true

// Define optional input
params.modules['star_genome_generate'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf" 
params.modules['star_genome_generate'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"

// Module inclusions
include { star_genome_generate } from '../main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

// Channel for FASTA file(s) 
Channel
    .value(["$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6a.fa", "$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6b.fa"])
    .set { ch_testData_fasta }

// Run tests
workflow {
    // Run genome indexing
    log.info ("Run genome indexing...")
    star_genome_generate ( params.modules['star_genome_generate'], ch_testData_fasta )

    // Collect file names and view output for single-end read mapping
    star_genome_generate.out.genomeIndex.collect() | view
    star_genome_generate.out.chrNameFile.collect() | view
    star_genome_generate.out.report.collect() | view

    // Check count of output files from star_genome_generate
    assert_channel_count( star_genome_generate.out.genomeIndex, "genomeIndex", 1 )
    assert_channel_count( star_genome_generate.out.chrNameFile, "chrNameFile", 1 )
    assert_channel_count( star_genome_generate.out.report, "report", 1 )
}