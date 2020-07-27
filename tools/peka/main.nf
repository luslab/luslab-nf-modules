#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition 
process peka {
    publishDir "${params.outdir}/peka",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-peka:latest'

    input:
      tuple val(sample_id), path(peaks), path(xls), path(genome), path(genome_index), path(regions)

    output:
        tuple val(sample_id), path("results/*.{pdf,tsv}"), emit: results

    script:

    """
    #!/usr/bin/env python
    import importlib.util
    spec = importlib.util.spec_from_file_location("peka", "/home/src/kmers.py")
    pe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pe)

    pe.run("$peaks",
     "$xls",
     "$genome",
     "$genome_index", 
     "$regions",
     ${params.peka_window},
     ${params.peka_window_distal},
     ${params.peka_kmer_length},
     ${params.peka_top_n},
     ${params.peka_percentile},
     ${params.peka_clusters},
     ${params.peka_smoothing},
     ${params.peka_all_outputs},
     ${params.peka_regions_selection},
     ${params.peka_subsample},
     "${params.peka_repeats}")
    """
}