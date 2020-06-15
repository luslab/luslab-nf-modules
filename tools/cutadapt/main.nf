#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "cutadapt/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'cutadapt'

// Specify DSL2
nextflow.preview.dsl = 2

// TODO check version of cutadapt in host process --> CUTADAPT 2.6 (latest is 2.9)

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'cutadapt'

//Sample IDs
params.internal_sample_id_enabled = false

/*-------------------------------------------------> CUTADAPT PARAMETERS <-----------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
CUSTOM PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

//Insert custom arguments
params.internal_custom_args = ''

/*-----------------------------------------------------------------------------------------------------------------------------
ADAPTER SEQUENCES PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/
//Single adapter sequence
//Automatically removes 3' adapters if not specified otherwise
params.internal_adapter_sequence = 'AGATCGGAAGAGC'

//Multiple adapter sequences
//Type the name of the FASTA file that contains the adapter sequences (file format -> name.fasta)
params.internal_multiple_adapters = false
params.internal_multi_adapt_fasta = ''

//Multiple adapter trimming options
params.internal_3_trim_multiple = false
params.internal_5_trim_multiple = false
params.internal_3_or_5_trim_multiple = false

/*-----------------------------------------------------------------------------------------------------------------------------
PAIRED-END READS PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

//Activate paired-end mode
params.internal_paired_end_mode = false 


/*-----------------------------------------------------------------------------------------------------------------------------
REPORTING PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/
//Suppress all output except error messages
params.internal_quiet = false

//Provides the GC content (as percentage) of the reads
params.internal_gc_content = 0

//Detailed information about where adapters were found in each read are written to the given file
params.internal_info_file = ''

//Report type changed to a one-line summary
params.internal_min_report = false

/*-----------------------------------------------------------------------------------------------------------------------------
BASIC PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

//Trims low-quality ends from reads
params.internal_min_quality = 10

//Disallow insertions and deletions entirely
params.internal_no_indels = false

//Determines the maximum error rate for a specific adaptor
params.internal_max_error_rate = 0

//Changes the minimum overlap length for all parameters
params.internal_min_overlap = 0

/*Unconditionally removes bases from the beginning or end of each read. If the given length is positive,
the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end.*/
params.internal_cut = 0

// ONLY AVAILABLE IN VERSION 2.8 AND ON SINGLE-END DATA
//Cutadapt searches both the read and its reverse complement for adapters
//params.internal_rev_comp = true

/*-----------------------------------------------------------------------------------------------------------------------------
SINGLE ADAPTER TRIMMING PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

//Remove 3' adapters --> DEFAULT OPTION
params.internal_3_trim = true

//Remove 5' adapters
params.internal_5_trim = false

//Can remove either 3' or 5' adapters
params.internal_3_or_5_trim = false

//Disallow internal matches for a 3’ adapter
params.internal_non_intern_3_trim = false

//Disallow internal matches for a 5' adapter
params.internal_non_intern_5_trim = false

//Anchor a 3' adapter to the end of the read
params.internal_anchor_3_trim = false

//Anchor a 5' adapter to the end of the read
params.internal_anchor_5_trim = false


/*-----------------------------------------------------------------------------------------------------------------------------
FILTERING PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

/*Discard processed reads that are shorter than LENGTH
To avoid issues with zero-length sequences, specify at least "-m 1" */
params.internal_min_length = 16

//Discard processed reads that are longer than LENGTH
params.internal_max_length = 0

//Discard reads in which an adapter was found.
params.internal_discard_trimmed = false

//Discard reads in which no adapter was found. This has the same effect as specifying --untrimmed-output /dev/null.
params.internal_discard_untrimmed = false

//Instead of discarding the reads that are too short according to -m, write them to FILE (in FASTA/FASTQ format).
params.internal_too_short_output = ''

//Instead of discarding reads that are too long (according to -M), write them to FILE (in FASTA/FASTQ format).
params.internal_too_long_output = ''

//Write all reads without adapters to FILE (in FASTA/FASTQ format) instead of writing them to the regular output file.
params.internal_untrimmed_output = ''

//Discard reads with more than COUNT N bases. If COUNT_or_FRACTION is a number between 0 and 1, it is interpreted as a fraction of the read length
params.internal_max_n_bases = 0

//Discard reads with more than ERRORS expected errors. The number of expected errors is computed as described in Edgar et al. (2015), (Section 2.2).
params.internal_max_expected_errors = 0

//Discard reads that did not pass CASAVA filtering. Illumina’s CASAVA pipeline in version 1.8 adds an is_filtered header field to each read. 
//Specifying this option, the reads that did not pass filtering (these are the reads that have a Y for is_filtered) will be discarded. Reads for which the header cannot be recognized are kept.
params.internal_discard_casava = false

/*-----------------------------------------------------------------------------------------------------------------------------
PAIRED-END READS TRIMMING PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

//Enables trimming of paired-end readings 
params.internal_paired_end_readings = false

/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

// Trimming reusable component
process cutadapt {
    tag "${sample_id}"

    publishDir "cutadapt/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fq.gz"), emit: trimmedReads
        path "*.txt", emit: report

    shell:
    
    //Building cutadapt arguments
    if (params.internal_custom_args == null){
        cutadapt_args = ''
    }else {
        cutadapt_args = "$params.internal_custom_args "
    }
    

    //Report types if-statements
    if (params.internal_quiet){
        cutadapt_args += "--quiet "
    }
    if (0 < params.internal_gc_content << 100){
        cutadapt_args += "--gc-content=$params.internal_gc_content "
    }
    if (params.internal_info_file != ''){
        cutadapt_args += "--info-file $params.internal_info_file "
    }
    if (params.internal_min_report){
        cutadapt_args += "--report=minimal "
    }

    //Basic parameters if-statements 
    if (params.internal_min_quality > 0){
        cutadapt_args += "-q $params.internal_min_quality "
    }
    if (params.internal_no_indels){
        cutadapt_args += "--no-indels "
    }
    if (params.internal_max_error_rate > 0){
        cutadapt_args += "-e $params.internal_max_error_rate "
    }
    if (params.internal_min_overlap > 0){
        cutadapt_args += "-0 ${params.internal_min_overlap} "
    }
    if (params.internal_cut != 0){
        cutadapt_args += "-u $params.internal_cut "
    }

    //Determines if there a single or multiple adapters
    if (params.internal_multiple_adapters == false){
        //Adapter trimming if-statements
        if (params.internal_3_trim){
            cutadapt_args += "-a $params.internal_adapter_sequence "
        }
        if (params.internal_5_trim){
            cutadapt_args += "-g $params.internal_adapter_sequence "
        }
        if (params.internal_3_or_5_trim){
            cutadapt_args += "-b $params.internal_adapter_sequence "
        }
        if (params.internal_non_intern_3_trim){
        X = "X"
            cutadapt_args += "-a $params.internal_adapter_sequence$X "
        }
        if (params.internal_non_intern_5_trim){
            cutadapt_args += "-g X$params.internal_adapter_sequence "
        }
        if (params.internal_anchor_3_trim){
            cutadapt_args += "-a $params.internal_adapter_sequence\$ "
        }
        if (params.internal_anchor_5_trim){
            cutadapt_args += "-g ^$params.internal_adapter_sequence "
        }
    }else {
        if (params.internal_3_trim_multiple){
            cutadapt_args += "-a file:$params.internal_multi_adapt_fasta "
        }
        if (params.internal_5_trim_multiple){
            cutadapt_args += "-g file:$params.internal_multi_adapt_fasta "
        }
        if (params.internal_3_or_5_trim_multiple){
            cutadapt_args += "-b file:$params.internal_multi_adapt_fasta "
        }
    } 

    //Filtering params
    if (params.internal_min_length > 0){
        cutadapt_args += "-m $params.internal_min_length "
    }
    if (params.internal_max_length != 0){
        cutadapt_args += "-M $params.internal_max_length "
    }
    if (params.internal_discard_trimmed){
        cutadapt_args += "--discard-trimmed "
    }
    if (params.internal_discard_untrimmed){
        cutadapt_args += "--discard-untrimmed "
    }
    if (params.internal_too_short_output != ''){
        cutadapt_args += "--too-short-output $params.internal_too_short_output "
    }
    if (params.internal_too_long_output != ''){
        cutadapt_args += "--too-long-output $params.internal_too_long_output "
    }
    if (params.internal_untrimmed_output != ''){
        cutadapt_args += "--untrimmed-output $params.internal_untrimmed_output "
    }
    if (params.internal_max_n_bases != 0){
        cutadapt_args += "--max-n $params.internal_max_n_bases "
    }
    if (params.internal_max_expected_errors != 0){
        cutadapt_args += "--max-expected-errors $params.internal_max_expected_errors "
    }
    if (params.internal_discard_casava){
        cutadapt_args += "--discard-casava $params.internal_discard_casava "
    }

    //Outputs and inputs
    cutadapt_args += "-o ${sample_id}.trimmed.fq.gz "

    //Paired-end mode -> determining paired output + inputs (forward and reverse)
    if (params.internal_paired_end_mode){
        internal_default_paired_end_args += "-p ${sample_id}.trimmed.fq.gz "
    }
    
    """
    cutadapt $cutadapt_args $reads > ${sample_id}_${params.internal_process_name}.txt
    """
}
