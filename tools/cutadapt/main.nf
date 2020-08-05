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
        tuple val(meta), path("*.fq.gz"), emit: fastq
        path "*.log", emit: report

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        // Construct CL line
        readList = reads.collect{it.toString()}
        if (readList.size > 1){
            cutadapt_command = "cutadapt ${args} -o ${prefix}.1.fq.gz -p ${prefix}.2.fq.gz $reads > ${meta.sample_id}_cutadapt.log"
        } else {
            cutadapt_command = "cutadapt ${args} -o ${prefix}.fq.gz $reads > ${meta.sample_id}_cutadapt.log"
        }
        // Log
        if (params.verbose){
            println ("[MODULE] cutadapt command: " + cutadapt_command)
        }

        //SHELL
        """
        ${cutadapt_command}
        """
}