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

// Compresses file to output (assumes running base system on linux)"
process compress {
    container 'ubuntu:16.04'

    input:
      tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("*.*"), emit: file
        path "*.*", emit: fileNoMeta

    script:
    """
    FILE=$input_file
    cat $input_file | gzip -v - > "\${FILE}.gz"
    """
}

// Generic awk process
process awk {
    publishDir "${params.outdir}/${opts.publish_dir}",
      mode: "copy", 
      overwrite: true,
      saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    container 'ubuntu:16.04'

    input:
      val opts
      tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("${outfile_name}"), emit: file
        path "${outfile_name}", emit: file_no_meta

    script:
        outfile_name = "awk_${input_file}"
        if(opts.outfile_name) {
          outfile_name = opts.outfile_name
        }

    """
    awk ${opts.args} $input_file > ${outfile_name}
    """
}
