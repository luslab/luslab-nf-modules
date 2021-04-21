#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for BLAST...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

// Two invocations of makeblastdb for two sequence types.
// The parameters for these are specified in blast.config - because if the
// publish directory for both of these tools is the same, it will overwrite them.
include {blast_makeblastdb as blast_makeblastdb_nucl} from "../main.nf"
include {blast_makeblastdb as blast_makeblastdb_prot} from "../main.nf"

// Search programs
include {blast_blastn} from "../main.nf"
include {blast_blastp} from "../main.nf"
include {blast_blastx} from "../main.nf"
include {blast_tblastn} from "../main.nf"
include {blast_tblastx} from "../main.nf"

// Utility programs
include {blast_windowmasker} from "../main.nf"
include {blast_asn_to_tab} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
test_reference_nucl = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/lambda1000a/lambda_top10.fasta"],
]
test_reference_prot = [
    [[sample_id:"ref_proteome"], "$baseDir/../../../test_data/lambda1000a/lambda_prot.faa"],
]

test_query_nucl = [
    [[sample_id:"query_nucl"], "$baseDir/../../../test_data/lambda1000a/lambda_holin.fna"],
]
test_query_prot = [
    [[sample_id:"query_prot"], "$baseDir/../../../test_data/lambda1000a/lambda_holin.faa"],
]

// Define test data input channels
Channel
    .from(test_reference_nucl)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_ref_nucl}
Channel
    .from(test_reference_prot)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_ref_prot}

Channel
    .from(test_query_nucl)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_nucl}
Channel
    .from(test_query_prot)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_query_prot}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Work
    blast_makeblastdb_nucl(params.modules['blast_makeblastdb_nucl'], ch_ref_nucl)
    blast_makeblastdb_prot(params.modules['blast_makeblastdb_prot'], ch_ref_prot)
    blast_blastn(params.modules['blast_blastn'], ch_ref_nucl, blast_makeblastdb_nucl.out.blast_db, ch_query_nucl)
    blast_blastp(params.modules['blast_blastp'], ch_ref_prot, blast_makeblastdb_prot.out.blast_db, ch_query_prot)
    blast_blastx(params.modules['blast_blastx'], ch_ref_prot, blast_makeblastdb_prot.out.blast_db, ch_query_nucl)
    blast_tblastn(params.modules['blast_tblastn'], ch_ref_nucl, blast_makeblastdb_nucl.out.blast_db, ch_query_prot)
    blast_tblastx(params.modules['blast_tblastx'], ch_ref_nucl, blast_makeblastdb_nucl.out.blast_db, ch_query_nucl)

    blast_asn_to_tab(params.modules['blast_asn_to_tab'], blast_blastn.out.asn, blast_makeblastdb_nucl.out.blast_db)
    blast_windowmasker(params.modules['blast_windowmasker'], ch_ref_nucl)

    // Print outputs
    blast_makeblastdb_nucl.out.blast_db | view
    blast_makeblastdb_prot.out.blast_db | view
    blast_blastn.out.asn | view
    blast_blastp.out.asn | view
    blast_blastx.out.asn | view
    blast_tblastn.out.asn | view
    blast_tblastx.out.asn | view
    blast_asn_to_tab.out.tab | view
    blast_windowmasker.out.fasta | view

    // Assert channel counts
    assert_channel_count(blast_makeblastdb_nucl.out.blast_db, "blast_db", 1)
    assert_channel_count(blast_makeblastdb_prot.out.blast_db, "blast_db", 1)
    assert_channel_count(blast_blastn.out.asn, "blastn", 1)
    assert_channel_count(blast_blastp.out.asn, "blastp", 1)
    assert_channel_count(blast_blastx.out.asn, "blastx", 1)
    assert_channel_count(blast_tblastn.out.asn, "tblastn", 1)
    assert_channel_count(blast_tblastx.out.asn, "tblastx", 1)
    assert_channel_count(blast_asn_to_tab.out.tab, "reformatted_blastn", 1)
    assert_channel_count(blast_windowmasker.out.wmstat, "windowmasker", 1)

}
