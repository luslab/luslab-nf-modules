#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Decompresses file to output (assumes running base system on linux)"
process decompress {
    input:
      tuple val(sample_id), path(input_file)

    output:
        tuple val(sample_id), path("*.*"), emit: file

    script: 
    """
    FILE=$input_file
    cat $input_file | gzip -dv - > "\${FILE%.*}"
    """
}