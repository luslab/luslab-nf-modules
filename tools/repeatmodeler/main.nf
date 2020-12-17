#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// This module uses the Dfam consortium's very handy container.
// See https://github.com/Dfam-consortium/TETools/ for how to add
// custom databases to the container.

// Process definition
process repeatmodeler_database {
    label "min_cores"
    label "min_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "dfam/tetools:1.2"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work --env "HOME=/work"'

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*{nhr,nin,nnd,nni,nog,nsq,translation}"), emit: repeatmodeler_db

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    repeatmodeler_database_command = "BuildDatabase $args -name ${fasta.simpleName} ${fasta}"

    if (params.verbose){
        println ("[MODULE] repeatmodeler_database_command command: " + repeatmodeler_database_command)
    }

    //SHELL
    """
    ${repeatmodeler_database_command}
    """
}

process repeatmodeler_model {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "dfam/tetools:1.2"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work --env "HOME=/work"'

    input:
        val opts
        tuple val(meta), path(fasta)
        tuple val(meta), path(repeatmodeler_db)

    output:
        tuple val(meta), path("*/consensi.fa"), emit: fasta
        tuple val(meta), path("*/families.stk"), emit: stockholm
        tuple val(meta), path("*"), emit: repeatmodeler

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    repeatmodeler_model_command = "RepeatModeler $args -srand ${opts.rng_seed} -database ${fasta.simpleName}"

    if (params.verbose){
        println ("[MODULE] repeatmodeler_model command: " + repeatmodeler_model_command)
    }

    //SHELL
    """
    ${repeatmodeler_model_command}
    """
}

process repeatclassifier {
    label "min_cores"
    label "min_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "dfam/tetools:1.2"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work --env "HOME=/work"'

    input:
        val opts
        tuple val(meta), path(fasta)
        tuple val(meta), path(stockholm)

    output:
        tuple val(meta), path("{*.classified,*.stk,*.fa}"), emit: report

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    repeatclassifier_command = "RepeatClassifier $args -engine ${opts.engine} -consensi ${fasta} -stockholm ${stockholm}"

    if (params.verbose){
        println ("[MODULE] repeatclassifier command: " + repeatclassifier_command)
    }

    //SHELL
    """
    ${repeatclassifier_command}
    """
}

process repeatmasker {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "dfam/tetools:1.2"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work --env "HOME=/work"'

    input:
        val opts
        tuple val(meta), path(query_fasta)
        tuple val(meta), path(lib_fasta)

    output:
        tuple val(meta), path("{*.out,*.tbl,*.log}"), emit: report
        tuple val(meta), path("*"), emit: repeatmasker

    script:

    args = ""
    if(opts.args) {
        ext_args = opts.args
        args += ext_args.trim()
    }

    // Optional arguments included by default here:
    //     -a      - write alignments to .align output file
    //     -html   - write optional HTML output file
    //     -gff    - write optional GFF output file
    //     -xsmall - soft-masking

    // Reminder: you can use a custom library with RepeatMasker by passing it
    // to the repeatmasker process as a meta/fasta channel.
    repeatmasker_command = "RepeatMasker $args -xsmall -a -html -gff -engine ${opts.engine} -pa ${opts.pa} -cutoff ${opts.cutoff} -lib ${lib_fasta} -dir ${query_fasta.simpleName}_masked ${query_fasta} > repeatmasker.log"

    if (params.verbose){
        println ("[MODULE] repeatmasker_command: " + repeatmasker_command)
    }

    //SHELL
    """
    ${repeatmasker_command}
    """
}
