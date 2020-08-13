#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*==============================
Test STAR genome indexing module
==============================*/

// Log
log.info ("Starting tests for STAR workflow: star_genome_generate -> star_align_reads...")

// Define main params
params.modules['star_genome_generate'].args = '--genomeSAindexNbases 9' 
params.modules['star_align_reads'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.verbose = true

// Define optional input
params.modules['star_genome_generate'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf"
params.modules['star_align_reads'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf" 
params.modules['star_genome_generate'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"
params.modules['star_align_reads'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"

// Module inclusions
include { star_genome_generate } from '../main.nf'
include { star_align_reads } from '../main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

// Channel for FASTA file(s) 
Channel
    .value(["$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6a.fa", "$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6b.fa"])
    .set { ch_testData_fasta }

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

// Run tests
workflow {
    // Run genome indexing and then read mapping
    log.info ("Run STAR workflow: star_genome_generate -> star_align_reads...")
    star_genome_generate ( params.modules['star_genome_generate'], ch_testData_fasta )
    star_align_reads ( params.modules['star_align_reads'], ch_testData_single_end, star_genome_generate.out.genomeIndex.collect() )   
    
    // Check count of output files from star_genome_generate
    assert_channel_count( star_genome_generate.out.genomeIndex, "genomeIndex", 1 )
    assert_channel_count( star_genome_generate.out.chrNameFile, "chrNameFile", 1 )
    assert_channel_count( star_genome_generate.out.report, "report", 1 )

    // Check count of output files from star_align_reads
    assert_channel_count( star_align_reads.out.samFiles, "samFiles", 2 )
    assert_channel_count( star_align_reads.out.bamFiles, "bamFiles", 2 )
    assert_channel_count( star_align_reads.out.sjFiles, "sjFiles", 2 )
    assert_channel_count( star_align_reads.out.chJunctions, "chJunctions", 0 )
    assert_channel_count( star_align_reads.out.readsPerGene, "readsPerGene", 0 )
    assert_channel_count( star_align_reads.out.finalLogFiles, "finalLogFiles", 2 )
    assert_channel_count( star_align_reads.out.outLogFiles, "outLogFiles", 2 )
    assert_channel_count( star_align_reads.out.progressLogFiles, "progressLogFiles", 2 )
    assert_channel_count( star_align_reads.out.report, "report", 2 )      
}
