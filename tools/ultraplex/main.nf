#!/usr/bin/env nextflow
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ULTRAPLEX {
    label "min_cores"
    label "min_memory"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${options.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else filename }
    
    container 'quay.io/biocontainers/ultraplex:1.1.4--py38h0213d0e_0'

    input:
        tuple val(meta), path(reads)
        val(barcode_file)

    output:
        tuple val(meta), path("*[!no_match].fastq.gz"), emit: fastq
        tuple val(meta), path("*no_match.fastq.gz"), optional: true, emit: no_match_fastq
        path "*.log", emit: report

    script:
        args = ""
        if(options.args && options.args != '') {
            ext_args = options.args
            args += ext_args.trim()
        }

        readList = reads.collect{it.toString()}
        if (readList.size > 1){
            ultraplex_command = "ultraplex --inputfastq ${readList[0]} --input_2 ${readList[1]} --barcodes $barcode_file --threads ${task.cpus} ${args}"
        } else {
            ultraplex_command = "ultraplex --inputfastq ${readList[0]} --barcodes $barcode_file --threads ${task.cpus} ${args}"
        }

        if (params.verbose){
            println ("[MODULE] ultraplex command: " + ultraplex_command)
        }

        """
        ${ultraplex_command}
        """
}
