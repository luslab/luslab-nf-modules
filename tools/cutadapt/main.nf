#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Trimming reusable component
process cutadapt {
    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }
    
    container 'luslab/nf-modules-cutadapt:latest'

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.trimmed.fq.gz"), emit: trimmedReads
        path "*.txt", emit: report

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        // Construct CL line
        cutadapt_command = "cutadapt ${args} -o ${meta.sample_id}.trimmed.fq.gz $reads > ${meta.sample_id}_cutadapt.txt"

        // Log
        if (params.verbose){
            println ("[MODULE] cutadapt command: " + cutadapt_command)
        }

        //SHELL
        readList = reads.collect{it.toString()}
        if (readList.size > 1){
            """
            cutadapt ${args} -o ${reads[0].simpleName}.trimmed.fq.gz -p ${reads[1].simpleName}.trimmed.fq.gz $reads > ${meta.sample_id}_cutadapt.txt
            """
        } else {
            """
            ${cutadapt_command}
            """
        }
}