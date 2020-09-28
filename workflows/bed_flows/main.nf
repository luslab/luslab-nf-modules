#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { bedtools_bamtobed } from '../../tools/bedtools/main.nf'
include { bedtools_genomecov } from '../../tools/bedtools/main.nf'
include { awk } from '../../tools/luslab_linux_tools/main.nf'
include { cut } from '../../tools/luslab_linux_tools/main.nf'
include { sort } from '../../tools/luslab_linux_tools/main.nf'
 
workflow paired_bam_to_bedgraph {
    take: bam
    take: genome
    take: meta
    main:

        // Define workflow parameters
        params.modules['bedtools_bamtobed'].args = '-bedpe'
        params.modules['awk'].args = '\'$1==$4 && $6-$2 < 1000 {print $0}\''
        params.modules['cut'].args = '-f 1,2,6'
        params.modules['sort'].args = '-k1,1 -k2,2n -k3,3n'
        params.modules['bedtools_genomecov'] = '-bg'


        // Convert BAM to BED
        bedtools_bamtobed( params.modules['bedtools_bamtobed'], meta, bam )

        // awk
        awk( params.modules['awk'], bedtools_bamtobed.out )

        // cut
        cut( params.modules['cut'], awk.out )

        // sort
        sort( params.modules['sort'], cut.out )

        // Get genome coverage in bedgraph format
        bedtools_genomecov( params.modules['bedtools_genomecov'], sort.out)


    emit:
        bedgraph = bedtools_genomecov.out.bed

}