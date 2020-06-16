#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "bedtools/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'bedtools_intersect'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'intersect'


// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

process bedtools_intersect {

    publishDir "bedtools/intersect/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input: 
        tuple val(sample_id), path(reads), path(regions_file)

    output: 
        tuple val(sample_id), path("${sample_id}.annotated.bed"), emit: annotatedBed

    shell:
    """
    bedtools intersect -a ${regions_file} -b $reads -wa -wb -s > ${sample_id}.annotated.bed
    """
}

