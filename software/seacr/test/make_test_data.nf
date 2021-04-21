#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include { paired_bam_to_bedgraph as hist_flow } from '../../../workflows/bed_flows/main.nf'
include { paired_bam_to_bedgraph as control_flow } from '../../../workflows/bed_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data paths

hist_data = [
    [[sample_id:"sample1_hist"], "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam.bai"],
    [[sample_id:"sample2_hist"], "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam.bai"]
]

control_data = [
    [[sample_id:"sample1_control"], "$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bam.bai"],
    [[sample_id:"sample2_control"], "$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bam.bai"]
]

genome_file = "$baseDir/../../../test_data/bed_flows/chr20.fa.fai"

// Define channels
Channel
    .from( hist_data )
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set { ch_hist_data }

Channel
    .from( control_data )
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set { ch_control_data }

Channel
    .value(file(genome_file))
    .set { ch_genome }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {

    // Run workflow for histone mod K27me3 data 

    hist_flow (ch_hist_data, ch_genome)

    // Run workflow for IgG control data 

    control_flow (ch_control_data, ch_genome)

}