#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for BUSCO...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {busco_genome} from "../main.nf"
include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_genome = [
	[[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/BUSCO_test_set_eukaryotes.fna"],
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
    busco_genome( params.modules["busco_genome"], ch_fasta )

    // Confirm the outputs of the above command
    busco_genome.out.report | view

    // Double check the channel count
		assert_channel_count( busco_genome.out.report, "busco", 1 )
}
