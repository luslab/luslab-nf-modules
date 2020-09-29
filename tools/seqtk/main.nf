#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Random subsample of FASTQ file
process seqtk_subsample {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-seqtk:base-1.0.0'
    
    input:
      val opts
      tuple val(meta), path(reads)

    output:
      tuple val(meta), path("*.fastq.gz"), emit: fastq
        
    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        //Seed parameter
        seed = "-s100"
        if(opts.seed && opts.seed != '') {
            seed = "-s" + opts.seed
        }

        //Number parameter
        number = 10000
        if(opts.base_count && opts.base_count != '') {
            number = opts.base_count
        }

        // Construct and log command
        seqtk_command = "seqtk sample $seed ${reads[0]} $number $args> ${prefix}.fastq"
        if (params.verbose){
            println ("[MODULE] seqtk subsample command: " + seqtk_command)
        }

        //Calculate the number of reads
        readList = reads.collect{it.toString()}
        if(readList.size == 1){
            """
            seqtk sample $seed ${reads[0]} $number $args> ${prefix}.fastq
            gzip ${prefix}.fastq
            """
        }
        else {
            """
            seqtk sample $seed ${reads[0]} $number $args> ${prefix}.r1.fastq
            seqtk sample $seed ${reads[1]} $number $args > ${prefix}.r2.fastq
            gzip ${prefix}.r1.fastq && gzip ${prefix}.r2.fastq
            """
        }
}

// Subset FASTA or FASTQ file with bed
process seqtk_subseq {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-seqtk:base-1.0.0'
    
    input:
      val opts
      path input
      path bed

    output:
        path "*${prefix}.*", emit: subset
        
    script:
        //Check extension
        ext = "fa"
        if ("$input" =~ /(fq.gz$)/) {
            ext = "fq"
        }

        prefix = opts.suffix ? "${input.simpleName}${opts.suffix}" : "${input.simpleName}"

        // Construct and log command
        seqtk_command = "seqtk subseq $input $bed > ${prefix}.${ext}"
        if (params.verbose){
            println ("[MODULE] seqtk subseq command: " + seqtk_command)
        }

        //SHELL
        """
        ${seqtk_command}
        """
}