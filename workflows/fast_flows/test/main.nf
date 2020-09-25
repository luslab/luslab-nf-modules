#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2
// Log
log.info ("Starting tests for fast_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {decompress} from '../../../tools/luslab_linux_tools/main.nf'
include {subset_genome} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .value([[:], "$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"])
    .set {ch_fasta}

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    decompress( ch_fasta )

    subset_genome( decompress.out.fileNoMeta, "21:40000000-40100000")
    subset_genome.out.fasta | view

    assert_channel_count( subset_genome.out.fasta, "fasta", 1)
}