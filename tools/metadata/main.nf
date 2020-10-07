#!/usr/bin/env nextflow

// Import groovy libs
import groovy.transform.Synchronized

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

workflow tenx_fastq_metadata {
    take: filePath
    main:
        Channel
            .fromPath( filePath )
            .splitCsv(header:true)
            .map { row -> processRow(row) }
            .map { row -> listFilesMultiDir(row, '.*.gz') }
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
            .map { row -> listFiles(row, '.*.gz') }
            .flatMap { row -> enumerateFastqDir(row) }
            .set { metadata }
    emit:
        metadata
}

@Synchronized
def listFiles(row, glob){
    file_array = []
    files = row[1].get(0).listFiles()
    for(def file:files){
        if(file.toString().matches(glob)){
            file_array.add(file)
        }
    }
    array = [row[0], [file_array]]
    return array
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

@Synchronized
def listFilesMultiDir(row, glob){
    file_array = []
    for(def dir:row[1]){
        for(def file:dir.listFiles()){
            if(file.toString().matches(glob)){
                file_array.add(file)
            }
        }
    }
    array = [row[0], file_array]
    return array
}

@Synchronized
def enumerateFastqDir(metadata){
    def fastqList = []
    def array = []
    if(metadata[0].strip2.isEmpty()){
        for (def fastq : metadata[1].flatten()){
            String s1 = fastq.getName().replaceAll(metadata[0].strip1, "")
            temp_meta = metadata[0].getClass().newInstance(metadata[0])
            temp_meta.remove("sample_id")
            temp_meta.put("sample_id", metadata[0].sample_id+"-"+s1)
            array.add([ temp_meta, [file(fastq, checkIfExists: true)]])
        }
    } else {
        fastqs = metadata[1].flatten().sort()

        for (int i = 0; i < fastqs.size(); i++){
            String s1 = fastqs.get(i).getName().replaceAll(metadata[0].strip1, "")
            temp_meta = metadata[0].getClass().newInstance(metadata[0])
            temp_meta.remove("sample_id")
            temp_meta.put("sample_id", metadata[0].sample_id+"-"+s1)
            array.add([ temp_meta, [file(fastqs.get(i), checkIfExists: true), file(fastqs.get(++i), checkIfExists: true)]])
        }
    }
    return array
}
