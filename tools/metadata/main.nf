#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

workflow fastq_metadata {
    take: filePath
    main:
        Channel
            .fromPath( filePath )
            .splitCsv(header:true)
            .map { row -> processRow(row) }
            .set { metadata }
    emit:
        metadata
}

def processRow(LinkedHashMap row) {
    def meta = [:]
    meta.sample_id = row.sample_id

    for (Map.Entry<String, ArrayList<String>> entry : row.entrySet()) {
        String key = entry.getKey();
        String value = entry.getValue();
    
        if(key != "sample_id" && key != "read1" && key != "read2") {
            meta.put(key, value)
        }
    }

    def array = []
    if (row.read2 == null) {
        array = [ meta, [ file(row.read1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true) ] ]
    }
    return array
}