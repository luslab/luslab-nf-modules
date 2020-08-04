#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*==============================
Test STAR genome indexing module
==============================*/

// Log
log.info ("Starting tests for STAR workflow: genomeGenerate -> alignReads...")

// Define main params
params.star_genomeGenerate_args = '--genomeSAindexNbases 9' 
params.star_alignReads_args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.verbose = true

// Define optional input
params.star_genomeGenerate_sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
params.star_alignReads_sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
//params.sjdbGTFfile = "$baseDir/input/hs_chr6_1Mbp/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf" //"gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
//"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
params.star_genomeGenerate_sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"
params.star_alignReads_sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"
params.genome_index = "${params.outdir}/star_genomeGenerate/genome_index"
//../../../test_data/star_index/hs_chr6_1Mb/2.7.5

// Module inclusions
include star_genomeGenerate from '../main.nf'
include star_alignReads from '../main.nf'

// Channel for FASTA file(s) 
Channel
    .fromPath("$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa")
    .set { ch_testData_fasta }

//     .combine( Channel.fromPath( params.genome_index ) )

// Single-end test reads
testMetaDataSingleEnd = [
  ['Sample1', "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep1.Unmapped.fq"],
  ['Sample2', "$baseDir/../../../test_data/fastq/prpf8_eif4a3_rep2.Unmapped.fq"]]

// Channel for single-end reads 
Channel
    .from( testMetaDataSingleEnd )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testData_single_end }

// Run tests
workflow {
    // Run genome indexing and then read mapping
    log.info ("Run STAR workflow: genomeGenerate -> alignReads...")
    star_genomeGenerate ( ch_testData_fasta )
    star_alignReads ( ch_testData_single_end
                          .combine( star_genomeGenerate.out.genomeIndex ) )
}