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

testDataGenome = [
	[[sample_id:"sample1"], "$baseDir/../../../test_data/fasta/Saccharomyces_cerevisiae_S288C-R64_chrI.fa"],
]

Channel
    .from(testDataGenome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run BUSCO on the test genome FASTA file
    busco_genome( params.modules["busco_genome"], ch_fasta )

    // Confirm the outputs of the above command
    busco_genome.out.busco | view
		busco_genome.out.report | view

    // Double check the channel count
		assert_channel_count( busco_genome.out.busco, "reads", 1 )
}
