#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bed_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include { paired_bam_to_bedgraph } from '../../../workflows/bed_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

bam_bai_test = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam", "$baseDir/../../../test_data/bed_flows/K27me3_1_to_chr20.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam", "$baseDir/../../../test_data/bed_flows/K27me3_2_to_chr20.bam.bai"]
]

genome_file = "$baseDir/../../test_data/bed_flows/chr20.genome"

// Define BAM channel
Channel
    .from(bam_bai_test)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set {ch_bam_bai}



/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    paired_bam_to_bedgraph ( ch_bam_bai, genome_file)

    paired_bam_to_bedgraph.out.bedgraph | view
}

