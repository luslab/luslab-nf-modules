#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process guppy_basecaller {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true
    
    if (params.num_gpus == 0){
    container "luslab/nf-modules-guppy:cpu"
    } else {
    container "luslab/nf-modules-guppy:gpu"
    }

    input:
        path(reads)

    output:
        tuple path("*.fastq"), path("*.log"), path("*.txt"), path("*.js"), emit: basecalledSeq
        path '*.txt', emit: summary
        
    script:

    // SHELL
    if (params.num_gpus == 0){
    """
    guppy_basecaller --input_path $reads --save_path . --flowcell ${params.guppy_flowcell} --kit ${params.guppy_kit}
    """
    } else {
    """
    guppy_basecaller --input_path $reads --save_path . --flowcell ${params.guppy_flowcell} --kit ${params.guppy_kit} -x cuda:all:100%
    """
    }
}

process guppy_qc {
    publishDir "${params.outdir}/guppy",
        mode: "copy", overwrite: true

    container "luslab/nf-modules-guppy:qc"

    input:
        path(summary) 

    output: 
        path '*.html', emit: report

    script:
    """
    pycoQC --summary_file $summary --html_outfile qc_report.html
    """
}