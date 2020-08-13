#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

//Process definition
process homer_annotatePeaks {
    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }
    
    container "quay.io/biocontainers/homer:4.11--pl526h2bce143_2"    

    input:
        val opts
        tuple val(meta), path(peaks)
        path(fasta)
        path(gtf)

    output:
        tuple val(meta), path("${prefix}"), emit: annotatedPeaks
    
    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

    // Construct CL line
    annotatePeaks_command = "annotatePeaks.pl ${args} ${peaks} ${fasta} -gtf ${gtf} > ${prefix}"

    // Log
    if (params.verbose){
        println ("[MODULE] homer/annotatePeaks command: " + annotatePeaks_command)
    }

    //SHELL
    """
    ${annotatePeaks_command}
    """
}

