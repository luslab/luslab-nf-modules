#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {r} from '../main.nf' addParams(options: modules['r'])

// Define test data channel
Channel.value(file("$baseDir/../../../test_data/r/test.csv"))
       .set {ch_test}

workflow {
    r (ch_test)
}