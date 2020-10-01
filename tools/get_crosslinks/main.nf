#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process getcrosslinks {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-get_crosslinks:base-1.0.0'

    input:
      val(opts)
      tuple val(meta), path(bam), path(bai)
      path(fai)

    output:
      tuple val(meta), path ("${prefix}.bed.gz"), emit: crosslinkBed

    script:

      prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

      //println ("[MODULE] getcrosslinks container: " + "${workflow.container}")

      //SHELL
      """
      bedtools bamtobed -i ${bam[0]} > dedupe.bed
      bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
      bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
      bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
      cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${prefix}.bed.gz
      """
}