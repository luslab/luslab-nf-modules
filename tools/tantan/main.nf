#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl = 2

// Process definition
process tantan {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "luslab/nf-modules-tantan:base-1.0.0"
    //container 'quay.io/biocontainers/tantan:13--he1b5a44_2'

    input:
        val opts
        tuple val(meta), path(genome)

    output:
        tuple val(meta), path("*.tantan.txt"), emit: tantanRepeats

    script:
	//Build the command line options
    command = "tantan -w${opts.max_period} -f4 ${genome} > ${meta.sample_id}.tantan.txt"
	//SHELL
    """
    ${command}
    """
}

process tantan_to_GFF3 {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'ubuntu:16.04'

    input:
        val opts
        tuple val(meta), path(tantan)

    output:
        tuple val(meta), path("*.tantan.gff3"), emit: tantanRepeats

    shell:
    '''
    awk '
        BEGIN {
            OFS="\t"
            print "##gff-version 3"
        }
        
        {
            print $1, "tantan", "tandem_repeat", $2 + 1, $3, int($4*$5), ".", ".", "Name="$6
        }
    ' !{tantan} > !{meta.sample_id}.tantan.gff3
    '''
}
