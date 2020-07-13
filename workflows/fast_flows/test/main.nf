#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for fast_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.genome = "$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz"
params.region = "21:40000000-40100000"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include decompress_noid from '../../../tools/luslab_file_tools/main.nf'
include subset_genome from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .from(params.genome)
    .set {ch_fasta}

//------------------------------------------------------------------------------------

// Run workflow
workflow {
    // Decompress fasta file
    decompress_noid( ch_fasta )

    // Call subset
    subset_genome( decompress_noid.out.file, params.region )

    // View output
    subset_genome.out.fastaSubset | view
}