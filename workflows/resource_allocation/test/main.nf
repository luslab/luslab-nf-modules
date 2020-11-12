#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2
// Log
log.info ("Starting tests for resource_allocation...")

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {assert_resource_allocation_models} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------

workflow {
    // Run the test script on resource allocation models
    assert_resource_allocation_models()
}