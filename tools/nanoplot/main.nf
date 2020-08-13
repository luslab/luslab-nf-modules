#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Trimming reusable component
process nanoplot {
  publishDir "${params.outdir}/nanoplot",
    mode: "copy", overwrite: true

  container "luslab/nf-modules-nanoplot:latest"

  input:
		val opts
    tuple val(sample_id), path(reads)
  output:
    tuple val(sample_id), path("*{pdf,html,log}"), emit: nanoplotOutputs
    path "*.html", emit: report

  script:

  // Construct command line
  nanoplot_command = "NanoPlot -t ${task.cpus} --N50 -f pdf -p ${sample_id}. --fastq $reads"

  // SHELL
  """
  ${nanoplot_command}
  """

}
