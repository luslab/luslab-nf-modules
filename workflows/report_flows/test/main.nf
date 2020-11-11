#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for report_flows...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include { meta_report_annotate } from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

report_meta_test = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/report_flows/bowtie2_stats_exp.txt"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/report_flows/bowtie2_stats_spike.txt"]
]

meta_bam_bai = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam", "$baseDir/../../../test_data/atac-seq/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam", "$baseDir/../../../test_data/atac-seq/sample2.bam.bai"]
]

// Define report channel
Channel
    .from( report_meta_test )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_report_meta }

Channel
    .from( meta_bam_bai )
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] }
    .set { ch_bami_bai_meta }

Channel
    .value("$baseDir/../../../assets/awk_scripts/bt2_report_to_csv.awk")
    .set {ch_awk_file_script}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Parse bt2 report 
    meta_report_annotate ( ch_report_meta, ch_bami_bai_meta, ch_awk_file_script, params.modules )
}