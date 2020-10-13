#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for file_tools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.modules['awk'].args = /-F "\t" 'BEGIN{OFS="\t";} ($10 <-2000) || ($10 >1500) && ($8!~\/^exon.\/) {print $0}'/
params.modules['cut'].args = "-f 2,4,5,7"
params.modules['sort'].args = "-k 2"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {decompress; compress; awk; cut; sort} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

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

cut_data = [
    [[:], "$baseDir/../../../test_data/homer/chr1.gtf"]
]

Channel
    .from(cut_data)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true)]]}
    .set {ch_cut}


//------------------------------------------------------------------------------------

// Run workflow
workflow {
    decompress ( ch_input )
	compress (decompress.out.file)
    awk ( params.modules['awk'], ch_awk )
    cut ( params.modules['cut'], ch_cut )
    sort ( params.modules['sort'], cut.out.file )

    decompress.out.file | view
    awk.out.file | view
    cut.out.file | view
    sort.out.file | view

    assert_channel_count( decompress.out.file, "file", 1)
    assert_channel_count( awk.out.file, "awk_file", 1)
    assert_channel_count( cut.out.file, "cut_file", 1)
    assert_channel_count( sort.out.file, "sort_file", 1)
    assert_channel_count( awk.out.file, "file", 1 )
    assert_channel_count( compress.out.file, "file", 1 )
}
