#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "get_crosslinks_coverage/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'getcrosslinkscoverage'

// Specify DSL2
nextflow.preview.dsl = 2

// TODO check version of cutadapt in host process --> CUTADAPT 2.6 (latest is 2.9)

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'getcrosslinkscoverage'

//Prefix to define the output file 
params.internal_output_prefix = ''

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

process getcrosslinkscoverage {
    publishDir "get_crosslinks_coverage/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(bed)

    output:
      tuple val(sample_id), path("${bed.simpleName}.bedgraph.gz"), emit: bedGraph
      tuple val(sample_id), path("${bed.simpleName}.norm.bedgraph.gz"), emit: normBedGraph

    script:
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
