#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Random subsample of FASTQ file
process seqtk_subsample {
    publishDir "${params.outdir}/seqtk",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-seqtk:latest'
    
    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("*.sub.fastq.gz"), emit: sampledReads
        
    script:
        //Custom args
        args = ""
        if(params.seqtk_subsample_args && params.seqtk_subsample_args != '') {
            args += params.seqtk_subsample_args + " "
        }

        //Seed parameter
        seed = "-s100"
        if(params.seqtk_subsample_seed && params.seqtk_subsample_seed != '') {
            seed = "-s" + params.seqtk_subsample_seed
        }

        //Number parameter
        number = 10000
        if(params.seqtk_subsample_number && params.seqtk_subsample_number != '') {
            number = params.seqtk_subsample_number
        }

        // Construct and log command
        seqtk_command = "seqtk sample $seed ${reads[0]} $number $args> ${reads[0].simpleName}.sub.fastq"
        if (params.verbose){
            println ("[MODULE] seqtk subsample command: " + seqtk_command)
        }

        //Calculate the number of reads
        readList = reads.collect{it.toString()}
        if(readList.size == 1){
            """
            ${seqtk_command}
            gzip ${reads[0].simpleName}.sub.fastq
            """
        }
        else {
            """
            ${seqtk_command}
            seqtk sample $seed ${reads[1]} $number $args > ${reads[1].simpleName}.sub.fastq
            gzip ${reads[0].simpleName}.sub.fastq && gzip ${reads[1].simpleName}.sub.fastq
            """
        }
}

// Subset FASTA or FASTQ file with bed
process seqtk_subseq {
    publishDir "${params.outdir}/seqtk",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-seqtk:latest'
    
    input:
      tuple path(input), path(bed)

    output:
        path "*.sub.*", emit: subsetFile
        
    script:
        //Check extension
        ext = "fa"
        if ("$input" =~ /(fq.gz$)/) {
            ext = "fq"
        }

        // Construct and log command
        seqtk_command = "seqtk subseq $input $bed > ${input.simpleName}.sub.${ext}"
        if (params.verbose){
            println ("[MODULE] seqtk subseq command: " + seqtk_command)
        }

        //SHELL
        """
        ${seqtk_command}
        """
}