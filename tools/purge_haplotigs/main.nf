#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process purge_haplotigs_hist {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "quay.io/biocontainers/purge_haplotigs:1.1.1--0"
    container "luslab/nf-modules-purge_haplotigs:1.0.0"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work'

    input:
        val opts
        tuple val(meta), path(bam)
        tuple val(meta), path(fasta)

    output:
        // The output of purge_haplotigs is a directory that is used repeatedly
        // and output files are added to this directory during subsequent steps.
        tuple val(meta), path("*"), emit: purge_haplotigs_hist

    script:
    //Build the command line options
    purge_haplotigs_command = "purge_haplotigs readhist -threads ${task.cpus} -b ${bam} -g ${fasta}"

    //SHELL
    """
    ${purge_haplotigs_command}
    """
}

process purge_haplotigs_minima {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "quay.io/biocontainers/purge_haplotigs:1.1.1--0"
    container "luslab/nf-modules-purge_haplotigs:1.0.0"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work'

    input:
        val opts
        tuple val(meta), path(bam)
        tuple val(meta), path(fasta)
        tuple val(meta), path(purge_haplotigs)

    output:
        tuple val(meta), path("low_mid_high.csv"), emit: csv
        tuple val(meta), path("critical_values.csv"), emit: csv2
        tuple val(meta), path("purge_haplotigs_minima.log"), emit: report

    script:
    purge_haplotigs_minima = "purge_haplotigs_minima.py /work/tmp_purge_haplotigs/MISC/${bam}.histogram.csv"

    //SHELL
    """
    ${purge_haplotigs_minima}
    """
}

process purge_haplotigs_contigcov {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "quay.io/biocontainers/purge_haplotigs:1.1.1--0"
    container "luslab/nf-modules-purge_haplotigs:1.0.0"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work'

    input:
        val opts
        tuple val(meta), path(bam)
        tuple val(meta), path(fasta)
        tuple val(meta), path(purge_haplotigs_hist)
        val(annotated_meta)

    output:
        tuple val(meta), path("coverage_stats.csv"), emit: csv

    script:
    // If cutoffs are not set by the user, try to guess them
    // from the shape of the histogram with the python script.
    if (opts.cutoff_low == "" && opts.cutoff_mid == "" && opts.cutoff_high == "" ){
        purge_haplotigs_contigcov = "purge_haplotigs contigcov -i ${bam}.gencov -j ${opts.junk} -s ${opts.suspect} -l ${annotated_meta.low} -m ${annotated_meta.mid} -h ${annotated_meta.high} -o coverage_stats.csv"
    } else {
        purge_haplotigs_contigcov = "purge_haplotigs contigcov -i ${bam}.gencov -j ${opts.junk} -s ${opts.suspect} -l ${opts.cutoff_low} -m ${opts.cutoff_mid} -h ${opts.cutoff_high} -o coverage_stats.csv"
    }

    //SHELL
    """
    ${purge_haplotigs_contigcov}
    """
}

process purge_haplotigs_purge {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    //container "quay.io/biocontainers/purge_haplotigs:1.1.1--0"
    container "luslab/nf-modules-purge_haplotigs:1.0.0"
    containerOptions '-u \$(id -u):\$(id -g) -v "$PWD":/work'

    input:
        val opts
        tuple val(meta), path(bam)
        tuple val(meta), path(fasta)
        tuple val(meta), path(purge_haplotigs_hist)
        tuple val(meta), path(purge_haplotigs_coverage)

    output:
        tuple val(meta), path("*"), emit: purge_haplotigs_purge

    script:
    // If a .bed file containing repeat locations is provided,
    // this can be used to help the splitting. Otherwise, the
    // max_match_cov parameter can be used to specify a maximum
    // coverage cutoff, as regions exceeding this value will be
    // assumed to be repetiive.
    if (opts.repeat_bed == ""){
        purge_haplotigs_purge = "purge_haplotigs purge -threads ${task.cpus} -d -align_cov ${opts.align_cov} -max_match ${opts.max_match_cov} -g ${fasta} -c ${purge_haplotigs_coverage} -b ${bam}"
    } else {
        purge_haplotigs_purge = "purge_haplotigs purge -threads ${task.cpus} -d -align_cov ${opts.align_cov} -repeats ${opts.repeat_bed} -g ${fasta} -c ${purge_haplotigs_coverage} -b ${bam}"
    }

    //SHELL
    """
    ${purge_haplotigs_purge}
    """
}

workflow purge_haplotigs {
    take: tuple_meta_bam
    take: tuple_meta_fasta
    main:
        // Get coverage histogram
        purge_haplotigs_hist(params.modules['purge_haplotigs'], tuple_meta_bam, tuple_meta_fasta)
        // Estimate critical values from histogram
        purge_haplotigs_minima(params.modules["purge_haplotigs"], tuple_meta_bam, tuple_meta_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist)

        // Channel operations to add the calculated minima to the metadata
        ch_split_path_meta = ch_bam
            .map { row -> [row[0].sample_id, row[1..-1]].flatten() }
        purge_haplotigs_minima.out.csv.splitCsv(header:true)
            .map { row -> [ row[0].sample_id, row[0] << row[1] ] }
            .join ( ch_split_path_meta )
            .map { row -> row[1..-1] }
            .set { ch_annotated_meta }

        // Get contig coverage values
        purge_haplotigs_contigcov(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist, ch_annotated_meta.flatten())
        // Actually perform the purging
        purge_haplotigs_purge(params.modules["purge_haplotigs"], ch_bam, ch_fasta, purge_haplotigs_hist.out.purge_haplotigs_hist, purge_haplotigs_contigcov.out.csv)
    emit:
        purge_haplotigs_purge.out.purge_haplotigs_purge
}
