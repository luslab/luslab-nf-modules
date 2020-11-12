#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { bedtools_genomecov_scale_bam } from '../../tools/bedtools/main.nf'
include { sort } from '../../tools/luslab_linux_tools/main.nf'
 
workflow paired_bam_to_bedgraph {
    take: tuple_meta_bam
    take: scale_factor
    main:

        params.modules['sort'].args = '-k1,1 -k2,2n'
        params.modules['bedtools_genomecov_scale_bam'].args = '-pc'

        // Get genome coverage in bedgraph format
        bedtools_genomecov_scale_bam( params.modules['bedtools_genomecov_scale_bam'], tuple_meta_bam, scale_factor )

        // Sort bedgraph
        sort( params.modules['sort'], bedtools_genomecov_scale_bam.out.bedgraph )

    emit:
        bedgraph = sort.out.file

}