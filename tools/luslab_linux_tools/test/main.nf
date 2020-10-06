#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for file_tools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.modules['awk'].args = /-F "\t" 'BEGIN{OFS="\t";} ($10 <-2000) || ($10 >1500) && ($8!~\/^exon.\/) {print $0}'/

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {decompress; compress; awk} from '../main.nf'
include {assert_channel_count as assert_channel_count_compress; assert_channel_count as assert_channel_count_awk} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testData = [
    [[:], "$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"]
]

Channel
    .from(testData)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_input}

awkData = [
    [[:], "$baseDir/../../../test_data/luslab_linux_tools/test.txt"]
]

Channel
    .from(awkData)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_awk}

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    decompress ( ch_input )
	compress (decompress.out.file)
    assert_channel_count_compress( compress.out.file, "file", 1 )

    awk ( params.modules['awk'], ch_awk )
    assert_channel_count_awk( awk.out.file, "file", 1 )
}
