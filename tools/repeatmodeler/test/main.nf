#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log out
log.info ("Starting tests for repeatmodeler...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true
// RepeatModeler uses a parallelized search strategy with 4 cores per "pa".
params.modules["repeatmodeler_model"].pa = 1
// Search for LTRs
params.modules["repeatmodeler_model"].args = "-LTRStruct"
// RepeatMasker optional arguments
// Search speeds: -s (slow) -q (quick) -qq (rush job)
params.modules["repeatmasker"].args = "-qq -cutoff 200 -nolow"
// Reminder: you can use a custom library with RepeatMasker by passing it
// to the repeatmasker process as a meta/fasta channel.

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {repeatmodeler_database} from "../main.nf"
include {repeatmodeler_model} from "../main.nf"
include {repeatclassifier} from "../main.nf"
include {repeatmasker} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

test_data_fasta = [
    [[sample_id:"test-sample"], "$baseDir/../../../test_data/fasta/S_cerevisiae_chrXII.fa"],
]

Channel
    .from(test_data_fasta)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_fasta}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run minionqc on the test set
    repeatmodeler_database(params.modules['repeatmodeler_database'], ch_fasta)
    repeatmodeler_model(params.modules['repeatmodeler_model'], ch_fasta, repeatmodeler_database.out.repeatmodeler_db)
    repeatclassifier(params.modules['repeatclassifier'], repeatmodeler_model.out.fasta, repeatmodeler_model.out.stockholm)
    repeatmasker(params.modules['repeatmasker'], ch_fasta, repeatmodeler_model.out.fasta)

    // Collect file names and view output
    repeatmodeler_database.out.repeatmodeler_db | view
    repeatmodeler_model.out.fasta | view
    repeatmodeler_model.out.stockholm | view
    repeatmodeler_model.out.repeatmodeler | view
    repeatclassifier.out.report | view
    repeatmasker.out.report | view
    repeatmasker.out.repeatmasker | view

    assert_channel_count(repeatmodeler_database.out.repeatmodeler_db, "repeatmodeler_db", 1)
    assert_channel_count(repeatmodeler_model.out.fasta, "repeatmodeler_model consensus", 1)
    assert_channel_count(repeatmodeler_model.out.stockholm, "repeatmodeler_model families stockholm", 1)
    assert_channel_count(repeatmodeler_model.out.repeatmodeler, "repeatmodeler_model everything else", 1)
    assert_channel_count(repeatclassifier.out.report, "repeatclassifier report", 1)
    assert_channel_count(repeatmasker.out.report, "repeatmasker report", 1)
    assert_channel_count(repeatmasker.out.repeatmasker, "repeatmasker everything else", 1)
}
