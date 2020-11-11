#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for guppy...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {guppy_basecaller} from '../main.nf'
include {guppy_qc} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/lambda1000a"]
]

Channel
    .from(test_data)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fast5}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run guppy_basecaller
    guppy_basecaller ( params.modules['guppy_basecaller'], ch_fast5 )
    guppy_qc ( params.modules['guppy_qc'], guppy_basecaller.out.sequencing_summary )

    // Collect file names and view output
    guppy_basecaller.out.fastq | view
    guppy_basecaller.out.log | view
    guppy_basecaller.out.sequencing_summary | view
    guppy_basecaller.out.telemetry | view
    guppy_qc.out.report | view

    assert_channel_count( guppy_basecaller.out.fastq, "basecaller", 1)
    assert_channel_count( guppy_qc.out.report, "qc_report", 1)
}
