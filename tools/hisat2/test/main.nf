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

include {hisat2_build; hisat2_align} from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/


testDataSingleEnd= [
    [[sample_id:'Sample1'], "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/se_1_trimmed.fq.gz"],
    [[sample_id:'Sample1'], "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/se_2.trimmed.fq.gz"]
]

testDataPairedEnd= [
    [[sample_id:'Sample1'], "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/pe_1a.trimmed.fq.gz", "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/pe_1b.trimmed.fq.gz"],
    [[sample_id:'Sample2'], "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/pe_2a.trimmed.fq.gz", "/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/pe_2b.trimmed.fq.gz"]
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
    .from("/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/Gallus_gallus.sub.fa")
    .set {ch_genome}

Channel
    .from("/Users/alex/dev/repos/luslab-nf-modules/test_data/hisat2/galgal6_splice_sites_subset.txt")
    .set {ch_splice_sites}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // Run hisat2
    //hisat2 ( ch_fastq_single_end )
    hisat2_build ( params.modules['hisat2'], ch_genome )
    hisat2_align ( params.modules['hisat2'], ch_fastq_paired_end, hisat2_build.out.genome_index, ch_splice_sites )

    // Collect file names and view output
    hisat2_align.out.sam | view
}