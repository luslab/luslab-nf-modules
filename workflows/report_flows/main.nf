#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { awk_file } from '../../tools/luslab_linux_tools/main.nf'

workflow meta_report_annotate {
    take: tuple_report_meta
    take: tuple_path_meta
    take: awk_script
    take: module_params
    main:

        // turn meta data to value from map
        ch_split_path_meta = tuple_path_meta
            .map { row -> [row[0].sample_id, row[1..-1]].flatten() }

        // Produce csv to encorporate into metadata
        awk_file ( module_params['awk_file'], tuple_report_meta, awk_script )

        // Encorporate report data into annotated meta and original data path
        awk_file.out.file
            .splitCsv(header:true)
            .map { row -> [ row[0].sample_id, row[0] << row[1] ] }
            .join ( ch_split_path_meta )
            .map { row -> row[1..-1] }
            .set { ch_annotated_meta }
        ch_annotated_meta | view

    emit : annotated_input = ch_annotated_meta
}