#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "get_crosslinks/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'getcrosslinks'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'getcrosslinks'

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

process getcrosslinks {
    publishDir "get_crosslinks/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(bam), path (fai)

    output:
      tuple val(sample_id), path ("${bam.simpleName}.xl.bed.gz"), emit: crosslinkBed

    script:
    """
    bedtools bamtobed -i $bam > dedupe.bed
    bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
    bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
    cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${bam.simpleName}.xl.bed.gz
    """
}