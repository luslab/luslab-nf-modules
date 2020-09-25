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

include {busco_genome} from '../main.nf'

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
    busco_genome( params.modules['busco_genome'], ch_fasta )
    busco_genome.out.busco_genome_out | view
}
