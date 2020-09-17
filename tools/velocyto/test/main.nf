#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for velocyto...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true


/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {velocyto_run_smartseq2 as velocyto_run_smartseq2_a; velocyto_run_smartseq2 as velocyto_run_smartseq2_b; velocyto_run_10x} from '../main.nf'
include {assert_channel_count} from "$baseDir/../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

bam = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/velocyto/S1_chr1_test.bam"],
    [[sample_id:'Sample2'], "$baseDir/../../../test_data/velocyto/S2_chr1_test.bam"]
]

cellranger_out = [
    [[sample_id:'Sample1'], "$baseDir/../../../test_data/velocyto/cellrangerOut"]
]

//Define test data input channels
Channel
    .from(bam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_sc_bam}

Channel
    .from(bam)
    .map { row -> [ [sample_id:"all_cells"], file(row[1], checkIfExists: true) ] }
    .groupTuple(by: 0)
    .set {ch_grouped_bam}

Channel
    .from(cellranger_out)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_cellranger_out}

Channel
    .value(file("$baseDir/../../../test_data/velocyto/chr1.gtf"))
    .set {ch_gtf}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/
  
workflow {
    // // Run velocyto for each cell indivudally (one output file per cell)
    // velocyto_run_smartseq2_a ( params.modules['velocyto_run_smartseq2'], ch_sc_bam, ch_gtf )
    // // Run velocyto for all cells grouped (single output file) - also output hdf5 file using '-d 1' argument
    // velocyto_run_smartseq2_b ( params.modules['velocyto_run_smartseq2_hdf5'], ch_grouped_bam, ch_gtf )

    velocyto_run_10x ( params.modules['velocyto_run_10x'], ch_cellranger_out, ch_gtf )

    // // Collect file names and view output
    // velocyto_run_smartseq2_a.out.velocyto | view
    // velocyto_run_smartseq2_b.out.velocyto | view
    // // velocyto_run_10x.out| view

    // //Check count
    // assert_channel_count( velocyto_run_smartseq2_a.out.velocyto, "velocyto", 2)
    // assert_channel_count( velocyto_run_smartseq2_b.out.velocyto, "velocyto", 1)
}