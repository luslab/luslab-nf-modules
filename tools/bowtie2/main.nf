#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process def
process bowtie2_align {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        val opts
        tuple val(meta), path(reads)
        path index

    output:
        tuple val(meta), path("*.sam"), emit: sam
        path "*stats.txt", emit: report

    script:
        args = "-p ${task.cpus} --no-unal"
        files = ''

        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        readList = reads.collect{it.toString()}
        if(readList.size > 1){
            files = '-1 ' + reads[0] + ' -2 ' + reads[1]
        }
        else {
            files = '-U ' + reads[0]
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        command = "bowtie2 -x ${index[0].simpleName} $args $files 2>bowtie2_stats.txt > ${prefix}.sam"
        if (params.verbose){
            println ("[MODULE] bowtie2 command: " + command)
        }

        """
        $command
        cat bowtie2_stats.txt
        """
}

process bowtie2_build {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        val opts 
        path fasta

    output:
        path "*.bt2*", emit: bowtieIndex
        path "*.log", emit: report

    script:
        // Check main args string exists and strip whitespace
        args = ''
        if(opts.args && opts.arg != '') {
            ext_args = opts.arg
            args += " " + ext_args.trim()
        }

         // Construct CL line
        command = "bowtie2-build${args} --threads ${task.cpus} $fasta ${fasta.simpleName} > summary.log"

        // Log
        if (params.verbose){
            println ("[MODULE] bowtie2-build command: " + command)
        }

        """
        $command
        """
}