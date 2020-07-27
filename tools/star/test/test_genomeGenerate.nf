#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

/*==============================
Test STAR genome indexing module
==============================*/

// Log
log.info ("Starting tests for STAR genome indexing...")

// Define main params
params.star_genomeGenerate_args = '--genomeSAindexNbases 9' 
params.verbose = true

// Define optional input
params.sjdbGTFfile = "$baseDir/../../../../hs_chr20/raw_genome/Homo_sapiens.GRCh38.100.chr20.gtf" //"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf"
//params.sjdbGTFfile = "$baseDir/input/hs_chr6_1Mbp/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf" //"gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
//"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
params.sjdbFileChrStartEnd = "$baseDir/../../../../hs_chr20/raw_genome/Sample1.SJ.out.tab" //"$baseDir/input/raw_genome/Sample1.SJ.out.tab"

// Module inclusions
include star_genomeGenerate from '../main.nf'

// Channel for FASTA file(s) 
Channel
    .fromPath("$baseDir/../../../../hs_chr20/raw_genome/chr20.fa")
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