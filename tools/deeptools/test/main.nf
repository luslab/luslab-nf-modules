#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for deeptools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.modules['deeptools_bam_pe_fragment_size'].suffix = '_dt_fragment'

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {deeptools_bam_pe_fragment_size} from '../main.nf'
include {decompress} from '../../luslab_linux_tools/main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
bam_bai_test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bam_bai/sample1.bam", "$baseDir/../../../test_data/bam_bai/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bam_bai/sample2.bam", "$baseDir/../../../test_data/bam_bai/sample2.bam.bai"]
]

blacklist_gz = "$baseDir/../../../test_data/blacklists/example_blacklist.bed"

// Channels
Channel
    .from(bam_bai_test_data)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)  ] }
    .set {ch_test_bam_bai}

Channel
    .from(blacklist_gz)
    .set { ch_blacklist }

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
workflow {

    // Run deeptools_bam_pe_fragment_size
    deeptools_bam_pe_fragment_size( params.modules['deeptools_bam_pe_fragment_size'], ch_test_bam_bai, ch_blacklist.collect() )

    // Collect file names and view output
    deeptools_bam_pe_fragment_size.out.fragment_size_summary | view
    deeptools_bam_pe_fragment_size.out.report_meta | view
    deeptools_bam_pe_fragment_size.out.report_no_meta | view
    deeptools_bam_pe_fragment_size.out.fragment_stats_meta | view

    // Check count
    assert_channel_count( deeptools_bam_pe_fragment_size.out.fragment_size_summary, "fragment_size_summary", 2 )
    assert_channel_count( deeptools_bam_pe_fragment_size.out.report_meta, "report", 2 )
    assert_channel_count( deeptools_bam_pe_fragment_size.out.report_no_meta, "report_no_meta", 2 )
    assert_channel_count( deeptools_bam_pe_fragment_size.out.fragment_stats_meta, "fragment_stats_meta", 2 )
}

