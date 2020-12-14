#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process deeptools_bam_pe_fragment_size {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    container 'quay.io/biocontainers/deeptools:3.5.0--py_0'

    input:
        val opts
        tuple val(meta), path(bam), path(bai)
        path(blacklist)

    output:
        tuple val(meta), path("${prefix}_raw.csv"), path("${prefix}_summary.csv"), emit: fragment_size_summary
        tuple val(meta), path("${prefix}_summary.csv"), emit: fragment_stats_meta
        tuple val(meta), path("${prefix}_log.txt"), emit: report_meta
        path "${prefix}_log.txt", emit: report_no_meta
        path "${prefix}_raw.csv", emit: fragment_no_meta

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        command = "bamPEFragmentSize -b ${bam} --outRawFragmentLengths ${prefix}_raw.csv --table ${prefix}_summary.csv -bl ${blacklist} ${args} > ${prefix}_log.txt"

    //SHELL
    """
    ${command}
    """
}