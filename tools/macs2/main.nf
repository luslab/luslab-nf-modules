#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// https://github.com/macs3-project/MACS
process macs2_callpeaks {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container 'quay.io/biocontainers/macs2:2.2.7.1--py37h516909a_0'

    input:
        val opts
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peaks
        tuple val(meta), path("*.xls")                   , emit: xls

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"
        args += '--format ' + opts.format

        command = "macs2 callpeak $args --gsize ${opts.gsize} --name $prefix --treatment ${bam[0]}"

        if (params.verbose){
            println ("[MODULE] macs2 command: " + command)
        }

        """
        $command
        """
}