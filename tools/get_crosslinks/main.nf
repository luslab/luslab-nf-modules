#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process getcrosslinks {
    label "avg_cores"
    label "high_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    // bedtools=2.29.0,pigz=2.3.4
    container 'quay.io/biocontainers/mulled-v2-c8623b4f6522dddd48913bd12dcf405d1d4f0ce1:10e4c359b727e884f6e19ee978f89c44dbaca255-0'

    input:
      val(opts)
      tuple val(meta), path(bam), path(bai)
      path(fai)

    output:
      tuple val(meta), path ("${prefix}.bed.gz"), emit: crosslinkBed

    script:

      prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

      //SHELL
      """
      bedtools bamtobed -i ${bam[0]} > dedupe.bed
      bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
      bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
      bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
      cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${prefix}.bed.gz
      """
}