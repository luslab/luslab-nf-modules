#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*===============
Test STAR modules
================*/

// Module inclusions
include { star_genome_generate as index_genome_only } from '../main.nf'
include { star_genome_generate as index_genome } from '../main.nf'
include { star_align_reads as map_se } from '../main.nf'
include { star_align_reads as map_pe } from '../main.nf'
include { star_align_reads as align_reads } from '../main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

// Define main params for both star_genome_generate and star_align_reads
// Genome indexing
params.modules['index_genome_only'].args = '--genomeSAindexNbases 9'
params.modules['index_genome'].args = '--genomeSAindexNbases 9' 
// Read mapping
params.modules['map_se'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.modules['map_pe'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.modules['align_reads'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM'
params.verbose = true

// Define optional input files for star_genome_generate
// GTF
params.modules['index_genome_only'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf"
params.modules['index_genome'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf"
// Splice junctions
params.modules['index_genome_only'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"
params.modules['index_genome'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"

// Define optional input files for star_align_reads
// GTF
params.modules['map_se'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf" 
params.modules['map_pe'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf"
params.modules['align_reads'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6a_chr6b.gtf"
// Splice junctions
params.modules['map_se'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"
params.modules['map_pe'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"
params.modules['align_reads'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1_chr6a_chr6b.SJ.out.tab"

// Define test data

// Channel for two FASTA files
Channel
    .value(["$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6a.fa", "$baseDir/../../../test_data/fasta/GRCh38.primary_assembly.genome_chr6b.fa"])
    .set { ch_testData_2fastas }

// Channel for one FASTA file
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

// Paired-end test reads
testMetaDataPairedEnd = [
  [[sample_id:'Sample1'], "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_1.fq.bz2", "$baseDir/../../../test_data/fastq/ENCFF282NGP_chr6_3400000_3500000_1000reads_2.fq.bz2"]
]

// Channel for paired-end reads
Channel
    .from( testMetaDataPairedEnd )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)] ] }
    .set { ch_testData_paired_end }

// Channel for genome index
Channel
    .value(file("$baseDir/../../../test_data/star_index/hs_chr6_1Mb/2.7.5a"))
    .set { ch_test_index_file }

// Run tests
workflow {

    // Test genome generation
    log.info ("Test star_genome_generate module alone...")

    index_genome_only ( params.modules['index_genome_only'], ch_testData_fasta )

    // Check count of output files from index_genome_only (star_genome_generate)
    assert_channel_count( index_genome_only.out.genomeIndex, "genomeIndex", 1 )
    assert_channel_count( index_genome_only.out.chrNameFile, "chrNameFile", 1 )
    assert_channel_count( index_genome_only.out.report, "report", 1 )
}

// Test single-end read alignment
/*log.info ("Test star_align_reads module with single-end reads...")

params.modules['star_align_reads'].publis_dir = "star_genome_generate_single_end"

workflow {

}


log.info ("Test star_align_reads module with paired-end reads...")

log.info ("Test a workflow with star_genome_generate and star_align_reads...")*/



// Run tests
/*workflow {
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
}*/
