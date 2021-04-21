#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process def
process velocyto_run_smartseq2 {
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
    
    container "quay.io/biocontainers/velocyto.py:0.17.17--py37h3c125cd_1"

    input:
        val opts
        tuple val(meta), path(reads)
        path gtf

    output:
        tuple val(meta), path("*[loomhdf5]"), emit: velocyto

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        velocyto_command = "velocyto run-smartseq2 ${args} -o . ${reads} ${gtf}"
        if (params.verbose){
            println ("[MODULE] velocyto/run_smartseq2 command: " + velocyto_command)
        }

        """
        ${velocyto_command}

        if [ -d dump ]; then
            mv dump/* .
        fi
        """
}