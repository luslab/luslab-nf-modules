#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process last_db {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(reference_sequences)

    output:
        tuple val(meta), path("*"), emit: last_db

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        last_command = "lastdb $args -P ${task.cpus} ${meta.sample_id} ${reference_sequences}"

        if (params.verbose){
            println ("[MODULE] last_db command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_train {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta_db), path(last_db)
        tuple val(meta), path(query_sequences)

    output:
        tuple val(meta), path("*.par"), emit: par

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        last_command = "last-train $args -P ${task.cpus} ${meta_db.sample_id} ${query_sequences} > ${meta_db.sample_id}__${meta.sample_id}.par"

        if (params.verbose){
            println ("[MODULE] last_train command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_align {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta_db), path(last_db)
        tuple val(meta_par), path(last_train_par)
        tuple val(meta), path(query_sequences)

    output:
        tuple val(meta), path("*.maf"), emit: maf
        tuple val(meta), path("*.tab"), emit: tab

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = "${meta_db.sample_id}__${meta.sample_id}"

        last_command = "lastal $args -P ${task.cpus} -p ${last_train_par} ${meta_db.sample_id} ${query_sequences} > ${prefix}.maf"
        last_postmask_command = "last-postmask ${prefix}.maf | maf-convert -n tab > ${prefix}.tab"

        if (params.verbose){
            println ("[MODULE] last_align command: " + last_command)
            println ("[MODULE] LAST postmask command: " + last_postmask_command)
        }

        //SHELL
        """
        ${last_command}
        ${last_postmask_command}
        """
}

process last_filter_one_to_one {
    label "min_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(unfiltered_maf)

    output:
        tuple val(meta), path("*.filtered.maf"), emit: maf
        tuple val(meta), path("*.tab"), emit: tab

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_filter_command = "last-split ${unfiltered_maf} | maf-swap | last-split | maf-swap > ${unfiltered_maf.simpleName}.filtered.maf"
        last_postmask_command = "last-postmask ${unfiltered_maf.simpleName}.filtered.maf | maf-convert -n tab > ${unfiltered_maf.simpleName}.filtered.tab"

        if (params.verbose){
            println ("[MODULE] LAST one-to-one filter command: " + last_filter_command)
            println ("[MODULE] LAST postmask command: " + last_postmask_command)
        }

        //SHELL
        """
        ${last_filter_command}
        ${last_postmask_command}
        """
}

process last_filter_one_to_many {
    label "min_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(unfiltered_maf)

    output:
        tuple val(meta), path("*.one-to-many.maf"), emit: maf

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        last_filter_one_to_many_command = "last-split ${unfiltered_maf} > ${unfiltered_maf.simpleName}.one-to-many.maf"

        if (params.verbose){
            println ("[MODULE] LAST one-to-many filter command: " + last_filter_one_to_many_command)
        }

        //SHELL
        """
        ${last_filter_one_to_many_command}
        """
}

process last_convert_maf {
    label "min_cores"
    label "low_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(maf)

    output:
        tuple val(meta), path("*.${opts.suffix}")

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        last_command = "maf-convert ${args} ${opts.suffix} ${maf} > ${maf.simpleName}.${opts.suffix}"

        if (params.verbose){
            println ("[MODULE] last_convert_maf command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_dotplot {
    label "min_cores"
    label "low_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1186--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(tab)

    output:
        tuple val(meta), path("*.${opts.suffix}"), emit: plot

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        last_command = "last-dotplot -x 2000 -y 2000 --sort1=1 --sort2=3 --strands2=1 --rot1=v --rot2=h ${tab} ${tab.simpleName}.${opts.suffix}"

        if (params.verbose){
            println ("[MODULE] last_dotplot command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}
