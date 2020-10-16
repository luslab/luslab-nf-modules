#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { awk_file } from '../../tools/luslab_linux_tools/main.nf'
include { processRow } from '../../tools/metadata/main.nf'

workflow meta_report_annotate {
    take: tuple_report_meta
    take: awk_script
    take: module_params
    main:

        // Produce csv to encorporate into metadata
        awk_file ( module_params['awk_file'], tuple_report_meta, awk_script )

        awk_file.out.file | view

        awk_file.out.file_no_meta
            .splitCsv(header:true)
            .view()

        // Encorporat

        // params.modules['awk'].args = awk_command

        // awk(params.modules['awk'], tuple_report_meta)

 
}