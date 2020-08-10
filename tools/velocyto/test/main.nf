#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for velocyto...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
params.samtools_view_region = '1'

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {samtools_index; samtools_view; samtools_faidx} from "$baseDir/../../samtools/main.nf"
include {velocyto_run_smartseq2} from '../main.nf'
include {assert_channel_count} from "$baseDir/../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

bam = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/velocyto/S1_chr1_test.bam"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/velocyto/S2_chr1_test.bam"]
]

//Define test data input channels

// for file in $(find /Volumes/lab-luscomben/home/users/thierya/analysis/ailin_scRNAseq/output/alignment_200423/ss11/3_bams/*.bam | head -2); do rsync -azP $file /Users/alex/dev/repos/luslab-nf-modules/test_data/velocyto ; done

//Single end
Channel
    .from(bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bam}

Channel
    .value(file("$baseDir/../../../test_data/velocyto/chr1.gtf"))
    .set {ch_gtf}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {

    velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_bam, ch_gtf )

    // Collect file names and view output
    velocyto_run_smartseq2.out.velocyto | view

    //Check count
    assert_channel_count( velocyto_run_smartseq2.out.velocyto, "x", 2)
}