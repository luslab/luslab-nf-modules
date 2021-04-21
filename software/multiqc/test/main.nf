#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for multiqc...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {multiqc} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .fromPath("$baseDir/../../../test_data/multiqc/fastqc/*.zip")
    .set {ch_fastqc}
    //.subscribe {log.info("$it")}

Channel
    .fromPath("$baseDir/../../../test_data/multiqc/cutadapt/*.log")
    .set {ch_cutadapt}
    //.subscribe {log.info("$it")}

Channel
    .fromPath("$baseDir/../../../test_data/multiqc/custom_config.yaml")
    .set {ch_config}
    //.subscribe {log.info("$it")}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    multiqc( params.modules['multiqc'], ch_config, ch_fastqc.mix(ch_cutadapt).collect())

    multiqc.out.report | view

    assert_channel_count( multiqc.out.report, "output", 1)
}
