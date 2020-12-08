#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for purge_haplotigs...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {purge_haplotigs_hist} from "../main.nf"
include {purge_haplotigs_minima} from "../main.nf"
include {purge_haplotigs_contigcov} from "../main.nf"
include {purge_haplotigs_purge} from "../main.nf"

// Note: using the reduced BAM file, the included python script will
// not be able to determine these values.
params.modules["purge_haplotigs"].cutoff_low = 14
params.modules["purge_haplotigs"].cutoff_mid = 72
params.modules["purge_haplotigs"].cutoff_high = 191
// "Junk" > 100 = nothing discarded
params.modules["purge_haplotigs"].junk = 101
// Alignment coverage cutoff of 60%
params.modules["purge_haplotigs"].align_cov = 60

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_bam = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/purge_haplotigs/cns_p_ctg.alignedChained.0.015.bam"],
]
test_data_fasta = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/purge_haplotigs/cns_p_ctg.fasta"],
]

//Define test data input channels
Channel
    .from(test_data_bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_bam}
Channel
    .from(test_data_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run the suggested purge_haplotigs workflow, step by step.
    purge_haplotigs_hist(params.modules["purge_haplotigs"], ch_bam, ch_fasta)
    purge_haplotigs_minima(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist)

    // Channel operations to add the calculated minima to the metadata
    ch_split_path_meta = ch_bam
        .map { row -> [row[0].sample_id, row[1..-1]].flatten() }
    purge_haplotigs_minima.out.csv.splitCsv(header:true)
        .map { row -> [ row[0].sample_id, row[0] << row[1] ] }
        .join ( ch_split_path_meta )
        .map { row -> row[1..-1] }
        .set { ch_annotated_meta }

    purge_haplotigs_contigcov(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist, ch_annotated_meta.flatten())
    purge_haplotigs_purge(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist, purge_haplotigs_contigcov.out.csv)

    // Collect file names and view output
    purge_haplotigs_hist.out.purge_haplotigs_hist | view
    purge_haplotigs_minima.out.csv | view
    purge_haplotigs_minima.out.csv2 | view
    purge_haplotigs_minima.out.report | view
    purge_haplotigs_contigcov.out.csv | view
    purge_haplotigs_purge.out.purge_haplotigs_purge | view

    // Verify channel counts
    //assert_channel_count(purge_haplotigs.out.purge_haplotigs, "whatever", 1)
}
