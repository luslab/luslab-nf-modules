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

workflow smartseq2_fastq_metadata {
    take: filePath
    main:
        Channel
            .fromPath( filePath )
            .splitCsv(header:true)
            .map { row -> processRow(row) }
            .flatMap { row -> enumerateFastqDir(row)}
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
    
        if(key != "sample_id" && key != "data1" && key != "data2") {
            meta.put(key, value)
        }
    }

    def array = []
    if (row.data2 == null) {
        array = [ meta, [ file(row.data1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.data1, checkIfExists: true), file(row.data2, checkIfExists: true) ] ]
    }
    return array
}

def enumerateFastqDir(metadata){
    def fastqList = []
    def array = []
    if(metadata[0].strip2.isEmpty()){
        for (def fastq : metadata[1].flatten()){
            String s1 = fastq.getName().replaceAll(metadata[0].strip1, "")
            metadata[0].put("cellID", s1)
            array.add([ metadata[0], [file(fastq, checkIfExists: true)]])
        }
    } else {
        fastqs = metadata[1].flatten().sort()

        for (int i = 0; i <fastqs.size(); i++ ){
            String s1 = fastqs.get(i).getName().replaceAll(metadata[0].strip1, "")
            metadata[0].put("cellID", s1)
            array.add([ metadata[0], [file(fastqs.get(i), checkIfExists: true), file(fastqs.get(++i), checkIfExists: true)]])
        }
    }
    return array
}