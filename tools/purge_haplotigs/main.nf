#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition

process purge_haplotigs_one_process {
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
        tuple val(meta), path("*"), emit: purge_haplotigs

    script:

    //Build the command line options
    purge_haplotigs_histo = "purge_haplotigs readhist -threads ${task.cpus} -b ${bam} -g ${fasta}"
    purge_haplotigs_minima = "python3 /purge_haplotigs_minima.py /work/tmp_purge_haplotigs/MISC/${bam}.histogram.csv"
    purge_haplotigs_contigcov = "purge_haplotigs -i /work/${bam}.gencov -l X -m Y -h Z"

    //SHELL
    """
    ${purge_haplotigs_histo}
    ${purge_haplotigs_minima}

    """
}



process purge_haplotigs {
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
        tuple val(meta), path("*"), emit: purge_haplotigs

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    purge_haplotigs_histo = "purge_haplotigs readhist $args -threads ${task.cpus}-b ${bam} -g ${fasta}"

    //SHELL
    """
    ${purge_haplotigs_histo}
    """
}


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
        tuple val(meta), path("*"), emit: purge_haplotigs

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    purge_haplotigs_command = "purge_haplotigs $args readhist -b ${bam} -g ${fasta}"

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
        //tuple val(meta), path("*"), emit: purge_haplotigs
        tuple

    script:
    purge_haplotigs_minima = "python3 /purge_haplotigs_minima.py /work/tmp_purge_haplotigs/MISC/${bam}.histogram.csv"

    //SHELL
    """
    ${purge_haplotigs_minima}
    """
}


process purge_haplotigs_minima_2 {
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
        tuple val(meta), path("*"), emit: purge_haplotigs

    script:
    // purge_haplotigs_minima = "python3 /purge_haplotigs_minima.py /work/tmp_purge_haplotigs/MISC/${bam}.histogram.csv"

    //SHELL
    """
    #!/usr/bin/env python3
    import sys
    from scipy.signal import argrelmin
    import numpy as np

    hist_dict = {}

    with open("/work/tmp_purge_haplotigs/MISC/!{bam}.histogram.csv") as i:
        for line in i:
            line = line.strip()
            # Keys are sequencing coverage values, values are counts of that coverage value
            hist_dict[line.split(',')[0]] = int(line.split(',')[1])

    local_minima = argrelmin(np.array(list(hist_dict.values())), order = 3)
    local_minima_x = [list(hist_dict.keys())[X] for X in local_minima[0]]
    print("%s\t%s\t%s" % (local_minima_x[0], local_minima_x[1], local_minima_x[2]))
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
        tuple val(meta), path(purge_haplotigs)

    output:
        tuple val(meta), path("*"), emit: purge_haplotigs

    script:
    //public class CSVParser { public String[] parse( String csvLine ) { def matcher = csvLine =~ /"([^"]*)"|(?<=,|^)([^,]*)(?:,|$)/ ; matcher.collect { it[1] } } }
    //CSVParser("/work/")

    purge_haplotigs_contigcov = "purge_haplotigs contigcov -i ${bam}.gencov -l "

    //SHELL
    """
    ${purge_haplotigs_minima}
    """
}
