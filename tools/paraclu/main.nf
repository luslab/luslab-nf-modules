#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process paraclu {
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    tag "${sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/paraclu:10--he513fc3_0'

    input:
        val opts
        tuple val(sample_id), path(crosslinks)

    output:
        tuple val(sample_id), path(output_name), emit: peaks

    script:
        output_name = "${sample_id}_paraclu_peaks_score_${opts.min_cluster_score}.bed.gz"

        // Adding 1 to the start position in the first awk converts to 1-based coordinates
        // Subtracting 1 from the start position in the second awk converts to 0-based BED format

        """
        gunzip -c $crosslinks | \
        awk '{OFS = "\t"}{print \$1, \$6, \$2+1, \$5}' | sort -k1,1 -k2,2 -k3,3n | \
        paraclu ${opts.min_cluster_score} - | \
        paraclu-cut -d ${opts.min_density_increase} -l ${opts.max_cluster_length} |
        awk '{OFS = "\t"}{print \$1, \$3-1, \$4, ".", \$6, \$2}' |
        sort -k1,1 -k2,2n | gzip > $output_name
        """
}