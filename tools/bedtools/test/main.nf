#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for bedtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.modules['bedtools_intersect_regions'].args = '-wa -wb -s'

params.modules['bedtools_intersect'].args = '-wa -wb'
params.modules['bedtools_intersect'].suffix = '+sample2.intersected.bed'

params.modules['bedtools_subtract'].args = '-A'
params.modules['bedtools_subtract'].suffix = '-sample2.subtracted.bed'

params.verbose = true
bedChannelExpected = 7



/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {bedtools_intersect_regions; bedtools_intersect; bedtools_subtract; bedtools_bamtobed; bedtools_genomecov; bedtools_genomecov_scale; bedtools_genomecov_bam} from '../main.nf'
include {assert_channel_count} from '../../../workflows/test_flows/main.nf'

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

bam_bai_test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam", "$baseDir/../../../test_data/atac-seq/sample1.bam.bai"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam", "$baseDir/../../../test_data/atac-seq/sample2.bam.bai"]
]

bam_test_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1.bam"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2.bam"]
]

genomecov_bed_data = [
    [[sample_id:"sample1"], "$baseDir/../../../test_data/atac-seq/sample1_bed.bed"],
    [[sample_id:"sample2"], "$baseDir/../../../test_data/atac-seq/sample2_bed.bed"]
]


genomecov_genome = "$baseDir/../../../test_data/atac-seq/genome.fa.fai"

// Define regions file input channel
Channel.value(file("$baseDir/../../../test_data/gtf/regions_GENCODE_v30.gtf.gz"))
       .set {ch_test_regions_file}

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_crosslinks}

// Define bam&bai test data input channel
Channel
    .from(bam_bai_test_data)
    .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)  ] }
    .set {ch_test_bam_bai}

// Define bed test data input channel
Channel
    .from(genomecov_bed_data)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_test_bed}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run bedtools_intersect_regions
    bedtools_intersect_regions (params.modules['bedtools_intersect_regions'], ch_test_crosslinks, ch_test_regions_file )
    // Run bedtools_intersect (this is to intersect two specific samples)
    bedtools_intersect (params.modules['bedtools_intersect'], ch_test_crosslinks.filter{it[0].sample_id == "sample1"}, ch_test_crosslinks.filter{it[0].sample_id == "sample2"}.map{it[1]} )
    // Run bedtools_subtract (this is to subtract one sample from another)
    bedtools_subtract (params.modules['bedtools_subtract'], ch_test_crosslinks.filter{it[0].sample_id == "sample1"}, ch_test_crosslinks.filter{it[0].sample_id == "sample2"}.map{it[1]} )
    // // Run bedtools_bamtobed
    bedtools_bamtobed (params.modules['bedtools_bamtobed'], ch_test_bam_bai)
    // Run bedtools_genomecov
    bedtools_genomecov (params.modules['bedtools_genomecov'], ch_test_bed, genomecov_genome)
    // Run bedtools_genomecov
    bedtools_genomecov_scale (params.modules['bedtools_genomecov_scale'], ch_test_bed_scale, genomecov_genome, 1)
    // Run bedtools_genomecov_bam
    bedtools_genomecov_bam (params.modules['bedtools_genomecov_bam'], ch_test_bam_bai)

    // Run bedtools_genomecov_bam

    // Collect file names and view output
    bedtools_intersect_regions.out.bed | view
    bedtools_intersect.out.bed | view
    bedtools_subtract.out.bed | view
    bedtools_bamtobed.out.bed | view
    bedtools_genomecov.out.bed | view
    bedtools_genomecov_scale.out.bed | view
    bedtools_genomecov_bam.out.bed | view

    //Check count
    assert_channel_count( bedtools_intersect_regions.out.bed, "bed", 6)
    assert_channel_count( bedtools_intersect.out.bed, "bed", 1)
    assert_channel_count( bedtools_subtract.out.bed, "bed", 1)
    assert_channel_count( bedtools_bamtobed.out.bed, "bed", 2)
    assert_channel_count( bedtools_genomecov.out.bed, "bed", 2)
    assert_channel_count( bedtools_genomecov_scale.out.bed, "bed", 2)
    assert_channel_count( bedtools_genomecov_bam.out.bed, "bed", 2)
}