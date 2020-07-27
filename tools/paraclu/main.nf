#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process paraclu {
    publishDir "${params.outdir}/paraclu",
        mode: "copy", overwrite: true
      
    container 'luslab/nf-modules-paraclu:latest'

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

    call(f'paraclu ${params.paraclu_min_value} "{paraclu_input}" | paraclu-cut.sh -l ${params.paraclu_max_cluster_length} -d ${params.paraclu_min_density_increase} > "{paraclu_output}"', shell=True)
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