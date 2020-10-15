#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process last_db {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        val(meta),  val(mode), path(fasta)

    output:
        tuple val(meta), path("*.lastdb.prj"), emit: last_db_prj
				tuple path("*.lastdb.bck"), path("*.lastdb.des"), path("*.lastdb.sds"), path("*.lastdb.ssp"), path("*.lastdb.suf"), path("*.lastdb.tis"), emit: last_db

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

				if( mode == "genome_vs_genome_near" )
            last_command = "lastdb $args -P ${tasks.cpus} -uNEAR -R01 ${fasta.simpleName} ${fasta}"
				else if ( mode == "genome_vs_genome_distant" )
            last_command = "lastdb $args -P ${tasks.cpus} -uMAM4 -R01 ${fasta.simpleName} ${fasta}"
				else if ( mode == "reads_vs_genome_simple_repeat_masking" )
            last_command = "lastdb $args -P ${tasks.cpus} -uNEAR -R01 ${fasta.simpleName} ${fasta}"
				else if ( mode == "reads_vs_genome_external_repeat_masking" )
            last_command = "lastdb $args -P ${tasks.cpus} -uNEAR -R11 -c ${fasta.simpleName} ${fasta}"
				else
            error "Invalid lastdb creation mode: ${mode}"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_train {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), val(mode), path(fastq)
        path(last_db_prj)

    output:
				tuple val(meta), path("*.par"), emit: last_train_par

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

				if( mode == "genome_vs_genome" )
            last_command = "last-train $args -P ${tasks.cpus} --revsym --matsym --gapsym -E0.05 -C2 ${last_db_prj.simpleName} ${fasta}"
				else if ( mode == "reads_vs_genome" )
            last_command = "last-train $args -P ${tasks.cpus} -Q 1 ${last_db_prj.simpleName} ${fastq} > ${last_db_prj.simpleName}.par"
				else
            error "Invalid lastdb training mode: ${mode}"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_align {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), val(mode), path(last_db_prj), path(last_train_par)
				path(query_fastq), optional: true
				path(query_fasta), optional: true

    output:
				tuple val(meta), path("*.maf"), emit: maf

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

				if( mode == "genome_vs_genome_near" )
            last_command = "lastal $args -P ${tasks.cpus} -m 50 -E0.05 -C2 -p ${last_train_par} ${last_db_prj.simpleName} ${query_fasta} | last-split -m1 > ${last_db_prj.simpleName}-${query_fasta.simpleName}.maf"
				else if( mode == "genome_vs_genome_distant" )
            last_command = "lastal $args -P ${tasks.cpus} -m 100 -E0.05 -C2 -p ${last_train_par} ${last_db_prj.simpleName} ${query_fasta} | last-split -m1 > ${query_fasta.simpleName}.maf"
				else if ( mode == "reads_vs_genome" )
            last_command = "lastal $args -P ${tasks.cpus} -p ${last_train_par} ${last_db_prj} ${query_fastq}"
				else
            error "Invalid lastdb training mode: ${mode}"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_filter_maf {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(unfiltered_maf)

    output:
				tuple val(meta), path("*.filtered.maf"), emit: maf
				path("*.tab"), emit: tab

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_rerarrange_command = "maf-swap ${unfiltered_maf} | last-split -m1 | maf-swap > ${unfiltered_maf.simpleName}.filtered.maf"
				last_postmask_command = "last-postmask ${unfiltered_maf.simpleName}.filtered.maf | maf-convert -n tab > ${unfiltered_maf.simpleName}.filtered.tab"

        if (params.verbose){
            println ("[MODULE] LAST one-to-one rearrangement command: " + last_rerarrange_command)
            println ("[MODULE] LAST postmask command: " + last_postmask_command)
        }

        //SHELL
        """
        ${last_rerarrange_command}
        ${last_postmask_command}
        """
}

process last_convert_maf_to_sam {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(maf)

    output:
				tuple val(meta), path("*.sam"), emit: sam

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = "maf-convert sam ${maf} > ${maf.simpleName}.sam"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}

process last_dotplot {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(tab)

    output:
        tuple val(meta), path("*.tiff"), emit: tiff

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_command = "last-dotplot -x 2000 -y 2000 --sort1=1 --sort2=3 --strands2=1 --rot1=v --rot2=h ${tab} ${tab.simpleName}.tiff"

        if (params.verbose){
            println ("[MODULE] last command: " + last_command)
        }

        //SHELL
        """
        ${last_command}
        """
}
