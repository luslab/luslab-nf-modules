#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process last_make_db {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

    input:
        val opts
        tuple val(meta), path(fasta_file)

    output:
        tuple val(meta), path("*.prj"), emit: last_db

    script:
        args = "-uNEAR -R11"
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = "lastdb $args -P${tasks.cpus} -c ${fasta_file.simpleName} ${fasta_file}"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_train_reads_on_genome {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

		input:
				val opts
				tuple val(meta), path(last_db_prj), path(fastq)

		output:
				tuple val(meta), path("*.par"), emit: last_train_params

		script:
				args = "-Q0"
				if(opts.args && opts.args != '') {
						ext_args = opts.args
						args += ext_args.trim()
				}

				prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

				last_command = "last-train $args -P ${tasks.cpus} ${last_db_prj.simpleName} ${fastq} > ${last_db_prj.simpleName}.par"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_align_reads_to_genome {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

    input:
        val opts
        tuple val(meta), path(last_train_params), path(last_db_prj), path(fastq)

    output:
        tuple val(meta), path("*.maf"), emit: last_alignment

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = "lastal -P ${task.cpus} -p ${last_train_params} ${last_db_prj.simpleName} ${fastq} | last-split > $prefix.maf "
        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_align_genome_to_genome {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

    input:
        val opts
        tuple val(meta), path(reads),
        path fasta_file

    output:
        tuple val(meta), path("*.maf"), emit: mappedReads

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = ""
        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_convert_maf_to_sam {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last"

    input:
        val opts
        tuple val(meta), path(reads),
        path fasta_file

    output:
        tuple val(meta), path("*.maf"), emit: mappedReads

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = ""
        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}
