#!/usr/bin/env nextflow

/*
========================================================================================
                         NF Modules Test Wrapper
========================================================================================
 #### Homepage / Documentation
 https://github.com/luslab/nf-modules
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for nf-modules...")

include { build_debug_param_summary; luslab_header} from './tools/luslab_util/main.nf'
include {assert_resource_allocation_models} from './workflows/resource_allocation/main.nf'

//------------------------------------------------------------------------------------

log.info luslab_header()
log.info build_debug_param_summary()

workflow {
    // Run the test script on resource allocation models
    assert_resource_allocation_models()
}