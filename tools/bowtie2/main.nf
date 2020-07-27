#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process def
process bowtie2_align {
    publishDir "${params.outdir}/bowtie2_align",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        tuple val(sample_id), path(reads), path(index)

    output:
        tuple val(sample_id), path("*.sam"), emit: alignedReads
        path "*stats.txt", emit: report

    script:
        args = "-p ${task.cpus} --no-unal"
        files = ''

        if(params.bowtie2_args && params.bowtie2_args != '') {
            ext_args = params.bowtie2_args
            args += ' ' + ext_args.trim()
        }

        readList = reads.collect{it.toString()}
        if(readList.size > 1){
            files = '-1 ' + reads[0] + ' -2 ' + reads[1]
        }
        else {
            files = '-U ' + reads[0]
        }

        command = "bowtie2 -x ${index[0].simpleName} $args $files 2>bowtie2_stats.txt > ${sample_id}.sam"
        if (params.verbose){
            println ("[MODULE] bowtie2 command: " + command)
        }

        """
        $command
        """
}

process bowtie2_build {
    publishDir "${params.outdir}/bowtie2_build",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        path fasta

    output:
        path "*.bt2*", emit: bowtieIndex
        path "*.log", emit: report

    script:
        // Check main args string exists and strip whitespace
        args = ''
        if(params.bowtie2_build_args && params.bowtie2_build_args != '') {
            ext_args = params.bowtie2_build_args
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