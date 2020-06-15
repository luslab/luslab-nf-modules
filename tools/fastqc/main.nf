#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "fastqc/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'fastqc'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'fastqc'

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

/*-------------------------------------------------> FASTQC PARAMETERS <-----------------------------------------------------*/

/*---------------------------------------------------------------------------------------------------------------------------*/

process fastqc {
  publishDir "fastqc/${params.internal_outdir}/${params.internal_process_name}",
    mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path ("*.{zip,html}"), emit: fastqcOutput
      path "*.{zip,html}", emit: report

    shell:
    reportname = "${params.internal_process_name}"
    if(params.internal_process_name != "fastqc") {
      reportname = "${params.internal_process_name}_fastqc"
    }

    """
    fastqc --threads ${task.cpus} $reads
    mv ${reads.simpleName}*.html ${sample_id}_${reportname}.html
    mv ${reads.simpleName}*.zip ${sample_id}_${reportname}.zip
    """
}