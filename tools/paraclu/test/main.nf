#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting paraclu module testing")

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

params.paraclu_min_value = 10
params.paraclu_max_cluster_length = 200
params.paraclu_min_density_increase = 2

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {paraclu} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Defining input channels
--------------------------------------------------------------------------------------*/

// Defining test data
testData = [
    ['Sample1', "$baseDir/../../../test_data/crosslinks/sample1.xl.bed.gz"],
    ['Sample2', "$baseDir/../../../test_data/crosslinks/sample2.xl.bed.gz"],
    ['Sample3', "$baseDir/../../../test_data/crosslinks/sample3.xl.bed.gz"],
    ['Sample4', "$baseDir/../../../test_data/crosslinks/sample4.xl.bed.gz"],
    ['Sample5', "$baseDir/../../../test_data/crosslinks/sample5.xl.bed.gz"],
    ['Sample6', "$baseDir/../../../test_data/crosslinks/sample6.xl.bed.gz"]
]

// Define test data input channels
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_crosslinks}

output_line_counts = [
    Sample1: 7,
    Sample2: 0,
    Sample3: 0,
    Sample4: 4,
    Sample5: 0,
    Sample6: 0
]

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

process count_lines {
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    tag "${sample_id}"

    container 'ubuntu:16.04'

    input:
        tuple val(sample_id), path(peaks)

    output:
        tuple val(sample_id), stdout, emit: line_count

    script:
    """
    echo -n "\$(gunzip -c $peaks | wc -l)"
    """
}

workflow {
    // Run paraclu
    paraclu ( ch_crosslinks )

    // Collect file names and view output
    paraclu.out.peaks | view

    count_lines( paraclu.out.peaks )

    count_lines.out.subscribe {
        if(output_line_counts[it[0]] != it[1].toInteger()) {
            throw new Exception("Sample " + it[0] + " is expected to have " + output_line_counts[it[0]] + " lines, but has " + it[1])
        }
    }
}