#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Define workflow to load up a csv file as a channel
workflow metadata {
    take: csv
    main:
        Channel
            .fromPath( csv )
            .splitCsv(header:true)
            .map { row -> [ row.sample_id, file(row.fastq, checkIfExists: true) ]  }
            .set { data }
    emit:
        data
}