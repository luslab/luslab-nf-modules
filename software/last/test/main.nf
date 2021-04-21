#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for LAST...")

/*-------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

// Different alignments require different alignment parameters.
// So, for each type of alignment, we can import the process
// as a placeholder, and fill in the "args" with the necessary
// parameters.

// Parameters for genome-genome alignment
// Two sets: for closely related sequences, and for distantly related
// sequences.
// The two DB commands do differ by parameter selection; see last.config.
include {last_db as last_db_genome_to_genome_near} from "../main.nf"
include {last_db as last_db_genome_to_genome_distant} from "../main.nf"

// Similarly the align commands differ by parameters
include {last_align as last_align_genome_to_genome_near} from "../main.nf"
include {last_align as last_align_genome_to_genome_distant} from "../main.nf"

// The train commands are actually the same, but they need to be named
// separately for Nextflow comptability.
include {last_train as last_train_genome_to_genome_near} from "../main.nf"
include {last_train as last_train_genome_to_genome_distant} from "../main.nf"

// As the train steps, these are the same processes imported to use
// different names  ¯\_(ツ)_/¯
include {last_filter_one_to_one as last_filter_one_to_one_near} from "../main.nf"
include {last_filter_one_to_one as last_filter_one_to_one_distant} from "../main.nf"
include {last_filter_one_to_many as last_filter_one_to_many_near} from "../main.nf"
include {last_filter_one_to_many as last_filter_one_to_many_distant} from "../main.nf"
include {last_dotplot as last_dotplot_near} from "../main.nf"
include {last_dotplot as last_dotplot_distant} from "../main.nf"

// Parameter sets for read-to-genome alignment:
include {last_db as last_db_reads_to_genome} from "../main.nf"
include {last_train as last_train_reads_to_genome} from "../main.nf"
include {last_align as last_align_reads_to_genome} from "../main.nf"
include {last_convert_maf as last_convert_maf_to_sam} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
test_reference_genome = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]
test_query_genome_near = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/last/Enterobacteria_phage_HK630.fasta"],
]
test_query_genome_distant = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/last/Shigella_sonnei_strain_6207_plasmid_unnamed1.fasta"],
]
test_query_reads = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fastq.gz"],
]

// Define test data input channels
Channel
    .from(test_reference_genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_ref_genome}
Channel
    .from(test_query_genome_near)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_genome_near}
Channel
    .from(test_query_genome_distant)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_genome_distant}
Channel
    .from(test_query_reads)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_reads}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Genome-genome alignment for near relatives
    // Make a reference genome database, train the genome-genome alignment parameters,
    // make the alignment, filter the alignment, then make a dot plot.
    last_db_genome_to_genome_near( params.modules['last_db_genome_to_genome_near'], ch_ref_genome )
    last_train_genome_to_genome_near( params.modules['last_train_genome_to_genome_near'], last_db_genome_to_genome_near.out.last_db, ch_query_genome_near )
    last_align_genome_to_genome_near( params.modules['last_align_genome_to_genome_near'], last_db_genome_to_genome_near.out.last_db, last_train_genome_to_genome_near.out.par, ch_query_genome_near)
    last_filter_one_to_many_near( params.modules['last_filter_one_to_many_near'], last_align_genome_to_genome_near.out.maf )
    last_filter_one_to_one_near( params.modules['last_filter_one_to_one_near'], last_align_genome_to_genome_near.out.maf )
    last_dotplot_near( params.modules['last_dotplot_near'], last_filter_one_to_one_near.out.maf )

    // Collect file names and view output
    last_db_genome_to_genome_distant( params.modules['last_db_genome_to_genome_distant'], ch_ref_genome )
    last_train_genome_to_genome_distant( params.modules['last_train_genome_to_genome_distant'], last_db_genome_to_genome_distant.out.last_db, ch_query_genome_distant )
    last_align_genome_to_genome_distant( params.modules['last_align_genome_to_genome_distant'], last_db_genome_to_genome_distant.out.last_db, last_train_genome_to_genome_distant.out.par, ch_query_genome_distant)
    last_filter_one_to_many_distant( params.modules['last_filter_one_to_many_distant'], last_align_genome_to_genome_distant.out.maf )
    last_filter_one_to_one_distant( params.modules['last_filter_one_to_one_distant'], last_align_genome_to_genome_distant.out.maf )
    last_dotplot_distant( params.modules['last_dotplot_distant'], last_filter_one_to_one_distant.out.maf )

    // Read-genome alignment for long reads
    // Make a reference genome database, train the read-genome alignment parameters,
    // make the alignment, then convert the maf to a sam file.
    last_db_reads_to_genome( params.modules['last_db_reads_to_genome'], ch_ref_genome )
    last_train_reads_to_genome( params.modules['last_train_reads_to_genome'], last_db_reads_to_genome.out.last_db, ch_query_reads )
    last_align_reads_to_genome( params.modules['last_align_reads_to_genome'], last_db_reads_to_genome.out.last_db, last_train_reads_to_genome.out.par, ch_query_reads )
    last_convert_maf_to_sam( params.modules['last_convert_maf_to_sam'], last_align_reads_to_genome.out.maf )

    // Double check the channel counts
    assert_channel_count( last_db_genome_to_genome_near.out.last_db, "last_db", 1 )
    assert_channel_count( last_db_genome_to_genome_distant.out.last_db, "last_db", 1 )
    assert_channel_count( last_db_reads_to_genome.out.last_db, "last_db", 1 )
    assert_channel_count( last_train_genome_to_genome_near.out.par, "last_train", 1 )
    assert_channel_count( last_train_genome_to_genome_distant.out.par, "last_train", 1 )
    assert_channel_count( last_train_reads_to_genome.out.par, "last_train", 1 )
    assert_channel_count( last_align_genome_to_genome_near.out.maf, "last_align", 1 )
    assert_channel_count( last_align_genome_to_genome_distant.out.maf, "last_align", 1 )
    assert_channel_count( last_align_reads_to_genome.out.maf, "last_align", 1 )
    assert_channel_count( last_filter_one_to_one_near.out.maf, "last_filter", 1 )
    assert_channel_count( last_filter_one_to_one_distant.out.maf, "last_filter", 1 )
    assert_channel_count( last_convert_maf_to_sam.out, "last_convert", 1 )
    assert_channel_count( last_dotplot_near.out.plot, "last_dotplot", 1 )
    assert_channel_count( last_dotplot_distant.out.plot, "last_dotplot", 1 )

}
