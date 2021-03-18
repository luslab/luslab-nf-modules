#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info ("Starting dsl2 test")

process activate_gtf {
    container "biocontainers/biocontainers:v1.2.0_cv1"

    input:
    path(gtf_file)

    output:
    path "genes.gtf", emit: gtf

    script:
    """
    echo "hello word"
    """

}

aws_gtf = "s3://ngi-igenomes/igenomes//Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
// aws_gtf = "/Users/westc/Nextflow/dev/repos/genes.gtf"

Channel
    .from(aws_gtf)
    .set {ch_gtf}

// def result = true

workflow {
    activate_gtf(ch_gtf)
    //activate_gtf.out.gtf | view

    activate_gtf.out.gtf.subscribe {
        Workflow.check_gtf_format( it, 30 )
    }

    log.info result.toString()
}
