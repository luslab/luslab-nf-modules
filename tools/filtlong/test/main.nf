#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for Infernal...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {filtlong} from "../main.nf"
params.modules["filtlong"].args = "--min_length 30000 --keep_percent 90"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_nanopore = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq.gz"],
]

Channel
    .from(test_data_nanopore)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_nanopore}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run BUSCO on the test genome FASTA file
    filtlong( params.modules["filtlong"], ch_nanopore )

    // Confirm the outputs of the above command
    filtlong.out.fastq | view

    // Double check the channel count
    assert_channel_count( filtlong.out.fastq, "filtered nanopore reads", 1 )
}
