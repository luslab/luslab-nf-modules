#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for samtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.samtools_index_args = ''
params.samtools_view_args = ''
params.samtools_view_region = 'chr21'

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

// include {samtools_index; samtools_view; samtools_faidx; samtools_sort} from '../main.nf'
// include {decompress_noid} from '../../../tools/luslab_file_tools/main.nf'
include {samtools_view2} from '../main.nf'
/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/bam_bai/sample1.bam", "$baseDir/../../../test_data/bam_bai/sample1.bam.bai"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/bam_bai/sample2.bam", "$baseDir/../../../test_data/bam_bai/sample2.bam.bai"]
]

Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_testDataView}

// Channel
//     .from( testData )
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set { ch_testDataIndex }

// Channel
//     .from( testData )
//     .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ]] }
//     .set { ch_testDataView }

// Channel
//     .from("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
//     .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run samtools index
    // samtools_index ( ch_testDataIndex )

    // // View output
    // samtools_index.out.baiFiles | view

    // // Run samtools index
    // samtools_view ( ch_testDataView )

    // samtools_sort ( params.modules['samtools_sort'], samtools_view.out )

    // // View output
    // samtools_sort.out.bam | view

    // //Test samtools faidx
    // decompress_noid( ch_fasta )
    // samtools_faidx( params.modules['samtools_faidx'], decompress_noid.out.file )
    // samtools_faidx.out.indexedFiles | view

    samtools_view2 ( params.modules['samtools_view'], ch_testDataView )
    samtools_view2.out | view
}
