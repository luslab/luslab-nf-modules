#!/usr/bin/env nextflow

// Include NfUtils
params.internal_classpath = "paraclu/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.internal_classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'paraclu'

// Specify DSL2
nextflow.preview.dsl = 2

// Define default nextflow internals
params.internal_outdir = 'results'
params.internal_process_name = 'paraclu'

// Local default params
params.internal_min_value = 10
params.internal_max_cluster_length = 200
params.internal_min_density_increase = 2

// Check if globals need to 
nfUtils.check_internal_overrides(module_name, params)

process paraclu {
    publishDir "paraclu/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(crosslinks)

    output:
      tuple val(sample_id), path("*_peaks.bed.gz"), emit: peaks

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from subprocess import call
    import sys
    import os

    file_in = "${crosslinks}"
    print(file_in)
    #file_in_name = file_in.name
    #print(file_in_name)

    if file_in.endswith('.gz'):
      file_out = file_in.replace('.bed', '_peaks.bed')
    else:
      file_out = file_in.replace('.bed', '_peaks.bed.gz')
    
    df_in = pd.read_csv(file_in,
                        names = ["chrom", "start", "end", "name", "score", "strand"],
                        header=None, sep='\t')

    df_out = df_in[['chrom', 'strand', 'start', 'score']]

    df_out.sort_values(['chrom', 'strand', 'start'], ascending=[True, True, True], inplace=True)

    paraclu_input = file_in + '.paraclu_input'
    paraclu_output = file_in + '.paraclu_output'

    df_out.to_csv(paraclu_input, sep='\t', header=None, index=None)

    call(f'paraclu ${params.internal_min_value} "{paraclu_input}" | paraclu-cut.sh -l ${params.internal_max_cluster_length} -d ${params.internal_min_density_increase} > "{paraclu_output}"', shell=True)
    df_in = pd.read_csv(paraclu_output,
                        names = ["sequence_name", "strand","start", "end", "number_of_positions",
                                "sum_of_data_values", "min_density", "max_density"],
                        header=None, sep='\t')
    df_in['fourth_column'] = '.'
    df_out = df_in[['sequence_name', 'start', 'end', 'fourth_column', 'sum_of_data_values', 'strand']]
    df_out.sort_values(['sequence_name','start', 'end', 'strand'],
                      ascending=[True, True, True, True], inplace=True)
    df_out.to_csv(file_out, sep='\t', header=None, index=None)
    call(f'rm "{paraclu_input}"', shell=True)
    call(f'rm  "{paraclu_output}"', shell=True)

    """
}