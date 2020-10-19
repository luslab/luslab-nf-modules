#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting tests for last...")

/*------------------------------------------------------------------------------------*/
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

// Parameter sets for genome-genome comparisons between close relatives:
include {last_db} from "../main.nf" as last_db_genome_to_genome_near
include {last_align} from "../main.nf" as last_align_genome_to_genome_near
params{
    modules {
        "last_db_genome_to_genome_near" {
            args             = "-uNEAR -R01"
            publish_dir      = "last_db_genome_to_genome_near"
            publish_results  = "all"
        }
        "last_align_genome_to_genome_near" {
            args             = "-m50 -E0.05 -C2"
            publish_dir      = "last_db_genome_to_genome_near"
            publish_results  = "all"
        }
    }
}
// Parameter sets for genome-genome comparisons between distant relatives:
include {last_db} from "../main.nf" as last_db_genome_to_genome_distant
params{
    modules {
        "last_db_genome_to_genome_distant" {
            args             = "-uMAM4 -R01"
            publish_dir      = "last_db_genome_to_genome_distant"
            publish_results  = "all"
        }
    }
}

include {last_train} from "../main.nf" as last_train_genome_to_genome
params{
    modules {
        "last_train_genome_to_genome" {
            args             = "--revsym --matsym --gapsym -E0.05 -C2"
            publish_dir      = "last_train_genome_to_genome"
            publish_results  = "all"
        }
    }
}

// Parameter sets for read-to-genome alignment:
include {last_db} from "../main.nf" as last_db_reads_to_genome_simple_masking
include {last_db} from "../main.nf" as last_db_reads_to_genome_external_masking
params{
    modules {
        "last_db_reads_to_genome_simple_masking" {
            args             = "-uNEAR -R01"
            publish_dir      = "last_db"
            publish_results  = "all"
        }
        "last_db_reads_to_genome_external_masking" {
            args             = "-uNEAR -R01"
            publish_dir      = "last_db"
            publish_results  = "all"
        }
    }
}

include {last_db} from "../main.nf" as last_db_genome_to_genome_distant
include {last_db} from "../main.nf" as last_db_reads_to_genome

// Import processes that do not change between
include {last_train} from "../main.nf"
include {last_align} from "../main.nf"
include {last_filter_maf} from "../main.nf"
include {last_convert_maf_to_sam} from "../main.nf"
include {last_dotplot} from "../main.nf"

include {assert_channel_count} from "../../../workflows/test_flows/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
test_alignment_files = [
    [[sample_id:"ref_genome"], "$baseDir/../../../test_data/last/E_coli_K-12.fna"],
]

// Define regions file input channel
Channel
    .fromPath("$baseDir/../../../test_data/fasta/homo-hg37-21.fa.gz")
    .set {ch_test_fasta}

// Define test data input channel
Channel
    .from(test_genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_ref_genome}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run last
    last_db.

    // Collect file names and view output


    // Double check the channel count
    assert_channel_count( last_db.out., "last_indices", 5 )
}
