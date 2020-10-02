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

include { paired_bam_to_bedgraph } from '../../../workflows/bed_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data paths

hist_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam.bai"]
]

control_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/IgG_1_to_chr20.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bam",  "$baseDir/../../../test_data/bed_flows/IgG_2_to_chr20.bam.bai"]
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
    set { ch_contorl_data }

Channel
    .value(file(genome_file))
    .set { ch_genome }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

// Run workflow for histone mod K27me3 data 

paired_bam_to_bedgraph (ch_hist_data, ch_genome)

// Run workflow for IgG control data 

paired_bam_to_bedgraph (ch_control_data, ch_genome)