#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { samtools_sort } from '../../tools/samtools/main.nf'
include { samtools_index } from '../../tools/samtools/main.nf'
 
workflow sort_index_bam {
    take: tuple_meta_bam
    main:
        // def modules {
        //     'samtools_sort' {
        //         args             = ""
        //         suffix           = "_sorted.bam"
        //         publish_dir      = "samtools_sort"
        //         publish_results  = "none"
        //     }
        //     'samtools_index' {
        //         args             = ""
        //         suffix           = ".bam.bai"
        //         publish_dir      = "samtools_index"
        //         publish_results  = "none"
        //     }
        // }

        // Define workflow parameters
        params.modules['samtools_sort'].publish_results = 'none'
        params.modules['samtools_index'].publish_results = 'none'

        // Sort bam
        samtools_sort( params.modules['samtools_sort'], tuple_meta_bam )

        // Index bam
        samtools_index( params.modules['samtools_index'], samtools_sort.out.bam )

    emit:
        bam_bai = samtools_index.out.bam
}