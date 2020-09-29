#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process guppy_basecaller {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container params.num_gpus == 0 ? "luslab/nf-modules-guppy:cpu-1.0.0" : "luslab/nf-modules-guppy:gpu-1.0.0"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.fastq.gz"), emit: fastq
        path "*.log", emit: log
        path "*.txt", emit: sequencing_summary
        path "*.js", emit: telemetry

    script:
        args = opts.args.trim()

        command = ""
        if (params.num_gpus == 0){
            command = "guppy_basecaller $args --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit} --num_callers ${opts.num_callers} --cpu_threads_per_caller ${task.cpus} --chunks_per_caller ${opts.chunks_per_caller} --chunk_size ${opts.chunk_size} --input_path $reads"
        } else {
            command = "guppy_basecaller $args --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit} --num_callers ${opts.num_callers} --cpu_threads_per_caller ${opts.cpu_threads_per_caller} --chunks_per_caller ${opts.chunks_per_caller} --chunk_size ${opts.chunk_size} --gpu_runners_per_device ${opts.gpu_runners_per_device} --input_path $reads -x cuda:all:100% "
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

    container "luslab/nf-modules-guppy:cpu-1.0.0"

    input:
        val opts
        path(sequencing_summary)

    output:
        path '*.html', emit: report

    script:
        prefix = opts.suffix ? "${opts.suffix}" : ""

        command = "pycoQC --summary_file $sequencing_summary --html_outfile ${prefix}qc_report.html"
        if (params.verbose){
            println ("[MODULE] guppy_qc command: " + command)
        }

        //SHELL
        """
        ${command}
        """
}
