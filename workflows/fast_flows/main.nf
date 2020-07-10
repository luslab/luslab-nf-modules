#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

include region2bed from '../tools/luslab_genome_tools/main.nf'

// Define workflow to subset and index a genome region fasta file
workflow subset_genome {
    take: filePath, region
    main:
        // Create channels from inputs
        Channel
            .fromPath( filePath, checkIfExists: true )
            .set { ch_genome }
        Channel
            .from( region )
            .set { ch_region }

        // Create bed from region
        region2bed( ch_region )
    //emit:
        //metadata
}