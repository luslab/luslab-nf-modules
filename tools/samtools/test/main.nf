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
include {decompress} from '../../../tools/luslab_linux_tools/main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testDataIndex = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/samtools/S1_chr1_test.bam"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/samtools/S2_chr1_test.bam"]
]

// testDataView = [
//     [[sample_id:'Sample1'], "$baseDir/../../../test_data/samtools/sample1.bam", "$baseDir/../../../test_data/samtools/sample1.bam.bai"],
//     [[sample_id:'Sample2'], "$baseDir/../../../test_data/samtools/sample2.bam", "$baseDir/../../../test_data/samtools/sample2.bam.bai"]
// ]

Channel
    .from( testDataIndex )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testDataIndex }

// Channel
//     .from( testDataView )
//     .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ]] }
//     .set { ch_testDataView }

Channel
    .value([[:], "$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"])
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run samtools index
    samtools_index ( params.modules['samtools_index'], ch_testDataIndex )
    samtools_index.out.bam | view
    assert_channel_count( samtools_index.out.bam, "bai", 2)

    // Run samtools view
    samtools_view ( params.modules['samtools_view'], ch_testDataIndex )
    samtools_view.out.bam | view
    assert_channel_count( samtools_view.out.bam, "bam", 2)

    // Run samtools sort
    samtools_sort ( params.modules['samtools_sort'], samtools_view.out.bam )
    samtools_sort.out.bam | view
    assert_channel_count( samtools_sort.out.bam, "bam", 2)

    //Test samtools faidx
    decompress( ch_fasta )
    samtools_faidx( params.modules['samtools_faidx'], decompress.out.file_no_meta )
    samtools_faidx.out.indexedFasta | view
    samtools_faidx.out.fai | view
    assert_channel_count( samtools_faidx.out.indexedFasta, "indexedFasta", 1)
    assert_channel_count( samtools_faidx.out.fai, "fai", 1)
}
