#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process blast_makeblastdb {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*"), emit: blast_db

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_makeblastdb_command = "makeblastdb $args -dbtype ${opts.dbtype} -in ${fasta} -out ${fasta.simpleName}"

        if (params.verbose){
            println ("[MODULE] blast_makeblastdb command: " + blast_makeblastdb_command)
        }

        //SHELL
        """
        ${blast_makeblastdb_command}
        """
}

process blast_blastn {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(ref_fasta)
        tuple val(meta), path(blast_db)
        tuple val(meta), path(query_fasta)

    output:
        tuple val(meta), path("*.asn"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_blastn_command = "blastn $args -evalue ${opts.evalue} -num_threads ${task.cpus} -query ${query_fasta} -db ${ref_fasta.simpleName} -out ${ref_fasta.simpleName}-${query_fasta.simpleName}.asn"

        if (params.verbose){
            println ("[MODULE] blast_blastn command: " + blast_blastn_command)
        }

        //SHELL
        """
        ${blast_blastn_command}
        """
}

process blast_blastp {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(ref_fasta)
        tuple val(meta), path(blast_db)
        tuple val(meta), path(query_fasta)

    output:
        tuple val(meta), path("*.asn"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_blastp_command = "blastp $args -evalue ${opts.evalue} -num_threads ${task.cpus} -query ${query_fasta} -db ${ref_fasta.simpleName} -out ${ref_fasta.simpleName}-${query_fasta.simpleName}.asn"

        if (params.verbose){
            println ("[MODULE] blast_blastp command: " + blast_blastp_command)
        }

        //SHELL
        """
        ${blast_blastp_command}
        """
}

process blast_blastx {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(ref_fasta)
        tuple val(meta), path(blast_db)
        tuple val(meta), path(query_fasta)

    output:
        tuple val(meta), path("*.asn"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_blastx_command = "blastx $args -evalue ${opts.evalue} -num_threads ${task.cpus} -query ${query_fasta} -db ${ref_fasta.simpleName} -out ${ref_fasta.simpleName}-${query_fasta.simpleName}.asn"

        if (params.verbose){
            println ("[MODULE] blast_blastx command: " + blast_blastx_command)
        }

        //SHELL
        """
        ${blast_blastx_command}
        """
}

process blast_tblastn {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(ref_fasta)
        tuple val(meta), path(blast_db)
        tuple val(meta), path(query_fasta)

    output:
        tuple val(meta), path("*.asn"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_tblastn_command = "tblastn $args -evalue ${opts.evalue} -num_threads ${task.cpus} -query ${query_fasta} -db ${ref_fasta.simpleName} -out ${ref_fasta.simpleName}-${query_fasta.simpleName}.asn"

        if (params.verbose){
            println ("[MODULE] blast_tblastn command: " + blast_tblastn_command)
        }

        //SHELL
        """
        ${blast_tblastn_command}
        """
}

process blast_tblastx {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(ref_fasta)
        tuple val(meta), path(blast_db)
        tuple val(meta), path(query_fasta)

    output:
        tuple val(meta), path("*.asn"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_tblastx_command = "tblastn $args -evalue ${opts.evalue} -num_threads ${task.cpus} -query ${query_fasta.simpleName} -db ${ref_fasta.simpleName} -out ${ref_fasta}-${query_fasta.simpleName}.asn"

        if (params.verbose){
            println ("[MODULE] blast_tblastx command: " + blast_tblastx_command)
        }

        //SHELL
        """
        ${blast_tblastx_command}
        """
}

process blast_asn_to_tab {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(blast_output)
        tuple val(meta), path(blast_db)

    output:
        tuple val(meta), path("*.tab"), emit: blast_output

    script:

        args = ""
        if(opts.args) {
            ext_args = opts.args
            args += ext_args.trim()
        }

        blast_asn_to_tab_command = "blast_formatter $args -max_target_seqs ${opts.max_target_seqs} -archive ${blast_output} -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore\" > ${blast_output.simpleName}.tab"

        if (params.verbose){
            println ("[MODULE] blast_asn_to_tab command: " + blast_asn_to_tab_command)
        }

        //SHELL
        """
        ${blast_asn_to_tab_command}
        """
}

process blast_windowmasker {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/blast:2.9.0--pl526he19e7b1_7"

    input:
        val opts
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.wmstat"), emit: wmstat
        tuple val(meta), path("*_masked.fasta"), emit: fasta

    script:

        blast_windowmasker_mkcounts_command = "windowmasker -mk_counts -in ${fasta} > ${fasta.simpleName}.wmstat"
        blast_windowmasker_mask_command = "windowmasker -ustat ${fasta.simpleName}.wmstat -in ${fasta} > ${fasta.simpleName}_masked.fasta"

        if (params.verbose){
            println ("[MODULE] blast_windowmasker_mkcounts command: " + blast_windowmasker_mkcounts_command)
            println ("[MODULE] blast_windowmasker_mask command: " + blast_windowmasker_mask_command)
        }

        //SHELL
        """
        ${blast_windowmasker_mkcounts_command}
        ${blast_windowmasker_mask_command}
        """
}
