#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Define workflow to load up a csv file as a channel
workflow simple_metadata {
    take: filePath
    main:
        Channel
            .fromPath( filePath )
            .splitCsv(header:true)
            .map { row -> [ row.sample_id, [ row.read1, row.read2 ]] }
            .set { ch_metadata }
    emit:
        ch_metadata
}