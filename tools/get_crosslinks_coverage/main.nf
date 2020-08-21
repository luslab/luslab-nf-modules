#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process getcrosslinkscoverage {
    publishDir "${params.outdir}/get_crosslinks_coverage",
        mode: "copy", overwrite: true
    
    container 'luslab/nf-modules-get_crosslinks_coverage'

    input:
        val(opt)
        tuple val(meta), path(bed)

    output:
        tuple val(meta), path("${bed.simpleName}.bedgraph.gz"), emit: bedGraph
        tuple val(meta), path("${bed.simpleName}.norm.bedgraph.gz"), emit: normBedGraph

    script:

    //SHELL
    """
    # Raw bedgraphs
    gunzip -c $bed | awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' | pigz > ${bed.simpleName}.bedgraph.gz
    
    # Normalised bedgraphs
    TOTAL=`gunzip -c $bed | awk 'BEGIN {total=0} {total=total+\$5} END {print total}'`
    echo \$TOTAL
    gunzip -c $bed | awk -v total=\$TOTAL '{printf "%s\\t%i\\t%i\\t%s\\t%f\\t%s\\n", \$1, \$2, \$3, \$4, 1000000*\$5/total, \$6}' | \
    awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' | \
    sort -k1,1 -k2,2n | pigz > ${bed.simpleName}.norm.bedgraph.gz
    """
}
