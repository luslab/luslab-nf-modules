#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bed_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include { paired_bam_to_bedgraph } from '../../../workflows/bed_flows/main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

bam_bai_test = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam", "$baseDir/../../../test_data/atac-seq/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam", "$baseDir/../../../test_data/atac-seq/sample2.bam.bai"]
]

genome_file = "$baseDir/../../../test_data/atac-seq/genome.fa.fai"

// Define BAM channel
Channel
    .from(bam_bai_test)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set {ch_bam_bai}

Channel
    .value(file(genome_file))
    .set {ch_genome}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run paired_bam_to_bedgraph
    paired_bam_to_bedgraph ( ch_bam_bai, ch_genome , 1)
    // View output
    paired_bam_to_bedgraph.out.bedgraph | view
    // Check count
    assert_channel_count( paired_bam_to_bedgraph.out.bedgraph, "bed", 2)
}

