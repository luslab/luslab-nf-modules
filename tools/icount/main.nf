#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "icount/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'icount'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'icount'

params.internal_output_prefix = '' //Prefix to define the output file 
params.internal_container = 'tomazc/icount:latest' // Set default container for running icount - override to specific custom container

/*-------------------------------------------------> iCOUNT PARAMETERS <-----------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
PEAK PARAMETERS
-------------------------------------------------------------------------------------------------------------------------------*/

params.internal_half_window = '3'
params.internal_fdr = '0.05'

/*-------------------------------------------------------------------------------------------------------------------------------*/

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

// Trimming reusable component
process icount {
    // Set tag to sample id
    tag "${sample_id}"

    publishDir "icount/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(bed), path(seg)

    output:
      tuple val(sample_id), path("${bed.simpleName}.xl.peaks.bed.gz"), emit: peaks
      tuple val(sample_id), path("${bed.simpleName}.scores.tsv"), emit: peak_scores
      tuple val(sample_id), path("${bed.simpleName}.xl.clusters.bed.gz"), emit: clusters

    shell:
    """
    iCount peaks $seg $bed ${bed.simpleName}.xl.peaks.bed.gz \
        --scores ${bed.simpleName}.scores.tsv \
        --half_window ${params.internal_half_window} \
        --fdr ${params.internal_fdr}

    zcat ${bed.simpleName}.xl.peaks.bed.gz | \
    bedtools merge -i stdin -s -d ${params.internal_half_window} -c 4,5,6 -o distinct,sum,distinct | \
    gzip > ${bed.simpleName}.xl.clusters.bed.gz
    """
}