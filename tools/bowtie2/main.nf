#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process def
process bowtie2_align {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    // bowtie2=2.4.1,samtools=1.9,pigz=2.3.4
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f480262c6feea34eb5a49c4fdfbb4986490fefbb-0'

    input:
        val opts
        tuple val(meta), path(reads)
        path index

    output:
        tuple val(meta), path("*.sam"), optional: true, emit: sam
        tuple val(meta), path("*.bam"), path("*.bai"), optional: true, emit: bam
        tuple val(meta), path("${prefix}${opts.unmapped_suffix}.1.fastq.gz"), path("${prefix}${opts.unmapped_suffix}.2.fastq.gz"), optional: true, emit: unmapped_fq_pe
        tuple val(meta), path("${prefix}${opts.unmapped_suffix}.fastq.gz"), optional: true, emit: unmapped_fq_s
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

        // If clause for creating unmapped filename if requested
        if(opts.unmapped_suffix && opts.unmapped_suffix != '') {
            if(readList.size > 1){
                args += ' --un-conc-gz ' + "${prefix}${opts.unmapped_suffix}" + '.fastq.gz'
            }
            else {
                args += ' --un-gz ' + "${prefix}${opts.unmapped_suffix}" + '.fastq.gz'
            }
        }

        // command = "bowtie2 -x ${index[0].simpleName} $args $files 2>${opts.summary_name}.txt > ${prefix}.sam"

        sort_command = "samtools sort -@ ${task.cpus} /dev/stdin > ${prefix}.bam"
        index_command = "samtools index -@ ${task.cpus} ${prefix}.bam"

        if ( opts.output_sam && opts.output_sam == true ) {
            command = "bowtie2 -x ${index[0].simpleName} $args $files 2>${opts.summary_name}.txt > ${prefix}.sam"
        }
        else {
            command = "bowtie2 -x ${index[0].simpleName} $args $files 2>${opts.summary_name}.txt | $sort_command && $index_command"
        }

        if (params.verbose){
            println ("[MODULE] bowtie2 command: " + command)
        }

        """
        $command
        cat ${opts.summary_name}.txt
        """
}

process bowtie2_build {
    tag "${fasta}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    // bowtie2=2.4.1,samtools=1.9,pigz=2.3.4
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f480262c6feea34eb5a49c4fdfbb4986490fefbb-0'

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