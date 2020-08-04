#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for hisat2...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {hisat2_build; hisat2_splice_sites; hisat2_splice_align as hisat2_splice_align_se; hisat2_splice_align as hisat2_splice_align_pe} from "$baseDir/../main.nf"
    
include {assert_channel_count as assert_channel_count_se; assert_channel_count as assert_channel_count_pe} from "$baseDir/../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/


testDataSingleEnd= [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/hisat2/se_1_trimmed.fq.gz"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/hisat2/se_2.trimmed.fq.gz"]
]

testDataPairedEnd= [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/hisat2/pe_1a.trimmed.fq.gz", "$baseDir/../../../test_data/hisat2/pe_1b.trimmed.fq.gz"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/hisat2/pe_2a.trimmed.fq.gz", "$baseDir/../../../test_data/hisat2/pe_2b.trimmed.fq.gz"]
]


//Define test data input channels

//Single end
Channel
    .from(testDataSingleEnd)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fastq_single_end}

// Paired-end
Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

Channel
    .from("$baseDir/../../../test_data/hisat2/Gallus_gallus.sub.fa")
    .set {ch_genome}

Channel
    .from("$baseDir/../../../test_data/hisat2/chr1.gtf")
    .set {ch_gtf}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {

    hisat2_build ( params.modules['hisat2'], ch_genome )
    hisat2_splice_sites ( params.modules['hisat2'], ch_gtf )

    // Run hisat2 for single end data
    hisat2_splice_align_se ( params.modules['hisat2'], ch_fastq_single_end, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )
    hisat2_splice_align_pe ( params.modules['hisat2'], ch_fastq_paired_end, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )

    // Collect file names and view output
    hisat2_splice_align_se.out.sam | view
    hisat2_splice_align_pe.out.sam | view

    //Check count
    assert_channel_count_se ( hisat2_splice_align_se.out.sam, "sam", 2)
    assert_channel_count_pe ( hisat2_splice_align_pe.out.sam, "sam", 2)
}