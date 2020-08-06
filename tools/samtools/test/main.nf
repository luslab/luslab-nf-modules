#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for samtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.verbose = true
params.modules['samtools_view'].args = "-b -h"

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/
include {samtools_index; samtools_view; samtools_faidx; samtools_sort} from '../main.nf'
include {decompress_noid} from '../../../tools/luslab_file_tools/main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/bam_bai/sample1.bam", "$baseDir/../../../test_data/bam_bai/sample1.bam.bai"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/bam_bai/sample2.bam", "$baseDir/../../../test_data/bam_bai/sample2.bam.bai"]
]

Channel
    .from( testData )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testDataIndex }

Channel
    .from( testData )
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ]] }
    .set { ch_testDataView }

Channel
    .from("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run samtools index
    samtools_index ( params.modules['samtools_index'], ch_testDataIndex )
    samtools_index.out.bai | view
    assert_channel_count( samtools_index.out.bai, "bai", 2)

    // Run samtools view
    samtools_view ( params.modules['samtools_view'], ch_testDataView )
    samtools_view.out.bam | view
    assert_channel_count( samtools_view.out.bam, "bam", 2)

    // Run samtools sort
    samtools_sort ( params.modules['samtools_sort'], samtools_view.out.bam )
    samtools_sort.out.bam | view
    assert_channel_count( samtools_sort.out.bam, "bam", 2)

    //Test samtools faidx
    decompress_noid( ch_fasta )
    samtools_faidx( params.modules['samtools_faidx'], decompress_noid.out.file )
    samtools_faidx.out.indexedFasta | view
    assert_channel_count( samtools_faidx.out.indexedFasta, "indexedFasta", 1)
}
