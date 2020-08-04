#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*==============================
Test STAR genome indexing module
==============================*/

// Log
log.info ("Starting tests for STAR genome indexing...")

// Define main params
params.star_genomeGenerate_args = '--genomeSAindexNbases 9' 
params.verbose = true

// Define optional input
params.star_genomeGenerate_sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
params.star_genomeGenerate_sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"

// Module inclusions
include { star_genomeGenerate } from '../main.nf'

// Channel for FASTA file(s) 
Channel
    .fromPath("$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa")
    .set { ch_testData_fasta }

// Run tests
workflow {
    // Run genome indexing
    log.info ("Run genome indexing...")
    star_genomeGenerate ( ch_testData_fasta )

    // Collect file names and view output for single-end read mapping
    star_genomeGenerate.out.genomeIndex.collect() | view
    star_genomeGenerate.out.chrNameFile.collect() | view
    star_genomeGenerate.out.report.collect() | view
}