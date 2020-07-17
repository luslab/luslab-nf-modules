#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

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

include {samtools_index; samtools_view; samtools_faidx} from '../main.nf'
include decompress_noid from '../../../tools/luslab_file_tools/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    ['Sample1', "$baseDir/../../../test_data/bam_bai/sample1.bam", "$baseDir/../../../test_data/bam_bai/sample1.bam.bai"],
    ['Sample2', "$baseDir/../../../test_data/bam_bai/sample2.bam", "$baseDir/../../../test_data/bam_bai/sample2.bam.bai"]
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
    samtools_index ( ch_testDataIndex )

    // View output
    samtools_index.out.baiFiles | view

    // Run samtools index
    samtools_view ( ch_testDataView )

    // View output
    samtools_view.out.bamFiles | view

    //Test samtools faidx
    decompress_noid( ch_fasta )
    samtools_faidx( decompress_noid.out.file )
    samtools_faidx.out.indexedFiles | view
}
