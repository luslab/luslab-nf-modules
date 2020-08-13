#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*======================
Test STAR mapping module
======================*/

// Log
log.info ("Starting tests for STAR mapping...")

// Define main params
params.modules['star_align_reads'].args = '--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM' 
params.verbose = true

// Define optional input
params.modules['star_align_reads'].sjdbGTFfile = "$baseDir/../../../test_data/gtf/gencode.v30.primary_assembly.annotation_chr6_34000000_35000000.gtf"
params.modules['star_align_reads'].sjdbFileChrStartEnd = "$baseDir/../../../test_data/star_splice_junctions/Sample1.SJ.out.tab"

// Module inclusions
include { star_align_reads as map_pe } from '../main.nf'
include { assert_channel_count } from '../../../workflows/test_flows/main.nf'

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
    map_pe ( params.modules['star_align_reads'], ch_testData_paired_end, ch_test_index_file )

    // Check count of output files from star_align_reads
    assert_channel_count( map_pe.out.samFiles, "samFiles", 1 )
    assert_channel_count( map_pe.out.bamFiles, "bamFiles", 1 )
    assert_channel_count( map_pe.out.sjFiles, "sjFiles", 1 )
    assert_channel_count( map_pe.out.chJunctions, "chJunctions", 0 )
    assert_channel_count( map_pe.out.readsPerGene, "readsPerGene", 0 )
    assert_channel_count( map_pe.out.finalLogFiles, "finalLogFiles", 1 )
    assert_channel_count( map_pe.out.outLogFiles, "outLogFiles", 1 )
    assert_channel_count( map_pe.out.progressLogFiles, "progressLogFiles", 1 )
    assert_channel_count( map_pe.out.report, "report", 1 )
}

