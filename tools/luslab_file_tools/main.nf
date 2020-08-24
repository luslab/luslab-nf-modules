#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Decompresses file to output (assumes running base system on linux)"
process decompress {
    input:
      tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("*.*"), emit: file

    script: 
    """
    FILE=$input_file
    cat $input_file | gzip -dv - > "\${FILE%.*}"
    """
}