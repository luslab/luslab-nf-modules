#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Decompresses file to output (assumes running base system on linux)"
process decompress {
    container 'ubuntu:16.04'

    input:
      tuple val(meta), path(input_file)
      
    output:
        tuple val(meta), path("*.*"), emit: file
        path "*.*", emit: fileNoMeta

    script: 
    """
    FILE=$input_file
    cat $input_file | gzip -dv - > "\${FILE%.*}"
    """
}

// Generic awk process

process awk {
    container 'ubuntu:16.04'

    input:
      val opts
      tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("${outfile_name}"), emit: file
        path "${outfile_name}", emit: fileNoMeta

    script:
        outfile_name = "awk_${input_file}"
        if(opts.outfile_name) {
          outfile_name = opts.outfile_name
        }

    """
    awk ${opts.args} $input_file > ${outfile_name}
    """
}

process cut {
      container 'ubuntu:16.04'

      input:
      val opts
      tuple val(meta), path(input_file)

      output:
      tuple val(meta), path("${outfile_name}"), emit: file
      path "${outfile_name}", emit: fileNoMeta

      script:
        outfile_name = "cut_${input_file}"
        if(opts.outfile_name) {
          outfile_name = opts.outfile_name
        }

      """
      cut ${opts.args} $input_file > ${outfile_name}
      """
}

process sort {
    container 'ubuntu:16.04'

    input:
      val opts
      tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("${outfile_name}"), emit: file
        path "${outfile_name}", emit: fileNoMeta

    script:
        outfile_name = "sort_${input_file}"
        if(opts.outfile_name) {
          outfile_name = opts.outfile_name
        }

    """
    sort ${opts.args} $input_file > ${outfile_name}
    """
}