#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process guppy_basecaller {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container params.num_gpus == 0 ? "luslab/nf-modules-guppy:cpu" : "luslab/nf-modules-guppy:gpu"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple path("*.fastq"), path("*.log"), path("*.txt"), path("*.js"), emit: basecalledSeq
        path '*.txt', emit: summary
        
    script:
        if (params.num_gpus == 0){
            """
            guppy_basecaller --input_path $reads --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit}
            """
        } else {
            """
            guppy_basecaller --input_path $reads --save_path . --flowcell ${opts.flowcell} --kit ${opts.kit} -x cuda:all:100%
            """
        }
}

process guppy_qc {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true

    container "luslab/nf-modules-guppy:cpu"

    input:
        path(summary) 

    output: 
        path '*.html', emit: report

    script:
    """
    pycoQC --summary_file $summary --html_outfile qc_report.html
    """
}