#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process minimap2 {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "luslab/nf-modules-minimap2:latest"
    container 'quay.io/biocontainers/minimap2:2.17--hed695b0_3'
    input:
        val opts
        tuple val(meta), path (reads)
        path fasta_file

    output:
        tuple val(meta), path ("*.paf"), emit: alignedReads

    script:

    args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        minimap2_command = "minimap2 $args ${fasta_file} $reads > ${prefix}.paf "
        
        if (params.verbose){
            println ("[MODULE] minimap2 command: " + minimap2_command)
        }

    // SHELL
    """
    ${minimap2_command}
    """
}