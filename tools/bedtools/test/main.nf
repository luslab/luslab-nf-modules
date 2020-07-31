#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bedtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.modules['bedtools_intersect'].args = '-wa -wb -s'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {bedtools_intersect} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
testData = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/bed/sample1.xl.bed.gz"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/bed/sample2.xl.bed.gz"],
    [[sample_id:"sample3"], "$baseDir/../../../test_data/bed/sample3.xl.bed.gz"],
    [[sample_id:"sample4"], "$baseDir/../../../test_data/bed/sample4.xl.bed.gz"],
    [[sample_id:"sample5"], "$baseDir/../../../test_data/bed/sample5.xl.bed.gz"],
    [[sample_id:"sample6"], "$baseDir/../../../test_data/bed/sample6.xl.bed.gz"]
]

// Define regions file input channel
Channel.value(file("$baseDir/../../../test_data/gtf/regions_GENCODE_v30.gtf.gz"))
       .set {ch_test_regions_file}

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_crosslinks}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run bedtools_intersect
    bedtools_intersect (params.modules['bedtools_intersect'], ch_test_crosslinks, ch_test_regions_file )

    // Collect file names and view output
    bedtools_intersect.out.bed | view
}