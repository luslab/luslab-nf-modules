#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for file_tools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.modules['awk'].args = "'{print NF}'"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {decompress; compress; awk} from '../main.nf'
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
    [[:], "$baseDir/../../../test_data/homer/chr1.gtf"]
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

    awk ( params.modules['awk'], ch_awk )

    decompress.out.file | view

    // assert_channel_count( decompress.out.file, "file", 1)
}
