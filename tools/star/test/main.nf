#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for STAR mapping...")

// Define main params
params.genome_index = "../hs_chr20/reduced_star_index/2.7.5a"
//params.genome_index = "tools/star/test/input/reduced_star_index/2.6.1c" //_SAnBases_7" // 2.7.5a 2.7.1a
params.star_map_args = '--outFilterMultimapNmax 20 --quantMode GeneCounts' //--quantMode TranscriptomeSAM GeneCounts' // '--quantMode GeneCounts'
params.verbose = true

// Define optional input
params.sjdbGTFfile = "$baseDir/../../../../hs_chr20/raw_genome/Homo_sapiens.GRCh38.100.chr20.gtf" //"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf"
//"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf"
//"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
params.sjdbFileChrStartEnd = '' //"$baseDir/../../../../hs_chr20/raw_genome/Sample1.SJ.out.tab" //"$baseDir/input/raw_genome/Sample1.SJ.out.tab"
params.varVCFfile = ''

// Module inclusions
include star_map as map_se from '../main.nf'
include star_map as map_pe from '../main.nf'

// Define input channels

// Single-end test reads
testMetaDataSingleEnd = [
  ['Sample1', "../hs_chr20/single_end/prpf8-hela-eif4a3-sirna-20190611-ju-2_trimmed_chr20_0_64444167.fq.gz"],
  ['Sample2', "../hs_chr20/single_end/prpf8-hela-eif4a3-sirna-20190611-ju-4_trimmed_chr20_0_64444167.fq.gz"]
]

//testMetaDataSingleEnd = [
//  ['Sample1', "$baseDir/input/single_end/prpf8_eif4a3_rep1.Unmapped.fq"],
//  ['Sample2', "$baseDir/input/single_end/prpf8_eif4a3_rep2.Unmapped.fq"]
//]

// Paired-end test reads
//testMetaDataPairedEnd = [
//  ['Sample1', "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/input/paired_end/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"]
//]

// Channel for single-end reads 
 Channel
    .from( testMetaDataSingleEnd )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData_single_end }

// Channel for paired-end reads
/*  Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .combine( Channel.fromPath( params.genome_index ) )
    .set { ch_testData_paired_end } */

// Run tests
workflow {
    // Run single-end read mapping
    log.info ("Run single-end read mapping...")
    map_se ( ch_testData_single_end )

    // Collect file names and view output for single-end read mapping
    map_se.out.bamFiles.collect() | view
    map_se.out.samFiles.collect() | view
    map_se.out.sjFiles.collect() | view
    map_se.out.finalLogFiles.collect() | view
    map_se.out.outLogFiles.collect() | view
    map_se.out.progressLogFiles.collect() | view 

    // Run paired-end read mapping
    /*log.info ("Run paired-end read mapping...")
    params.paired_end = true
    map_pe ( ch_testData_paired_end )

    // Collect file names and view output for paired-end read mapping
    map_pe.out.bamFiles.collect() | view
    map_pe.out.samFiles.collect() | view
    map_pe.out.sjFiles.collect() | view
    map_pe.out.finalLogFiles.collect() | view
    map_pe.out.outLogFiles.collect() | view
    map_pe.out.progressLogFiles.collect() | view */
    //map_pe.out.readsPerCount.collect() | view
}