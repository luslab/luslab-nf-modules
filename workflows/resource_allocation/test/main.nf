#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2
// Log
log.info ("Starting tests for resource_allocation...")

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include { build_debug_param_summary; luslab_header} from '../../../tools/luslab_util/main.nf'
include {assert_resource_allocation_models} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------

log.info luslab_header()
def summary = [:]
summary['Max CPUs'] = params.max_cpus
summary['Max memory'] = params.max_memory
summary['Max time'] = params.max_time
summary['Max GPUs'] = params.max_gpus
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

workflow {
    // Run the test script on resource allocation models
    assert_resource_allocation_models()
}