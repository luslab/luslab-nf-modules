#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

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
params.sjdbGTFfile = "$baseDir/../../../../hs_chr20/raw_genome/Homo_sapiens.GRCh38.100.chr20.gtf" //"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf"
//params.sjdbGTFfile = "$baseDir/input/hs_chr6_1Mbp/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000_first_gene.gtf" //"gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
//"$baseDir/input/raw_genome/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf" 
params.sjdbFileChrStartEnd = "$baseDir/../../../../hs_chr20/raw_genome/Sample1.SJ.out.tab" //"$baseDir/input/raw_genome/Sample1.SJ.out.tab"
params.genome_index = "${params.outdir}/star_genomeGenerate/genome_index"

// Module inclusions
include star_genomeGenerate from '../main.nf'
include star_alignReads from '../main.nf'

// Channel for FASTA file(s) 
Channel
    .fromPath("$baseDir/../../../../hs_chr20/raw_genome/chr20.fa")
    .set { ch_testData_fasta }

//     .combine( Channel.fromPath( params.genome_index ) )

// Single-end test reads
testMetaDataSingleEnd = [
  ['Sample1', "../hs_chr20/single_end/prpf8-hela-eif4a3-sirna-20190611-ju-2_trimmed_chr20_0_64444167.fq.gz"],
  ['Sample2', "../hs_chr20/single_end/prpf8-hela-eif4a3-sirna-20190611-ju-4_trimmed_chr20_0_64444167.fq.gz"]
]

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