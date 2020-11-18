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

include {infernal_cmscan} from "../main.nf"
params.modules["infernal_cmscan"].search_space = "1000000"
params.modules["infernal_cmscan"].db = "$baseDir/../../../test_data/infernal/RF00001.cm"
params.modules["infernal_cmscan"].clanin = "$baseDir/../../../test_data/infernal/RF00001.clanin"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_genome = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/S_cerevisiae_chrXII.fa"],
]

Channel
    .from(test_data_genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run BUSCO on the test genome FASTA file
    infernal_cmscan( params.modules["infernal_cmscan"], ch_fasta )

    // Confirm the outputs of the above command
    infernal_cmscan.out.tbl | view

    // Double check the channel count
    assert_channel_count( infernal_cmscan.out.tbl, "cmscan", 1 )
}
