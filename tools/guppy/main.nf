#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process guppy_basecaller {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container params.num_gpus == 0 ? "luslab/nf-modules-guppy:cpu" : "luslab/nf-modules-guppy:gpu"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.fastq.gz"), emit: fastq
        path "*.log", emit: log
        tuple path("*.txt"), path("*.js"), emit: report

    script:
        args = opts.args.trim()

        command = ""
        if (params.num_gpus == 0){
            command = "guppy_basecaller $args --input_path $reads --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit} --num-callers ${task.cpus} --cpu_threads_per_caller ${opts.threads_per_caller} --records_per_fastq 0"
        } else {
            command = "guppy_basecaller $args --input_path $reads --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit} --num-callers ${task.cpus} --records_per_fastq 0 -x cuda:all:100%"
        }

        if (params.verbose){
            println ("[MODULE] guppy command: " + command)
        }

        """
        $command
        """
}

process guppy_qc {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-guppy:cpu"

    input:
        val opts
        tuple path(summary), path(telemetry)

    output: 
        path '*.html', emit: report

    script:
        prefix = opts.suffix ? "${opts.suffix}" : ""

        command = "pycoQC --summary_file $summary --html_outfile ${prefix}qc_report.html"
        if (params.verbose){
            println ("[MODULE] guppy_qc command: " + command)
        }

        //SHELL
        """
        ${command}
        """
}
