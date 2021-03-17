#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Log
log.info ("Starting dsl2 test")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

process activate_gtf {
    input:
    path(gtf_file)

    output:
    path("genes.gtf"), emit: gtf

    script:
    """
    echo "hello word"
    """

}

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

aws_gtf = "s3://ngi-igenomes/igenomes//Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
// aws_gtf = "/Users/westc/Nextflow/dev/repos/genes.gtf"

// Define BAM channel
Channel
    .from(aws_gtf)
    .set {ch_gtf}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    activate_gtf(ch_gtf)

    activate_gtf.out.gtf | view
        
    gtf_check = check_gtf_by_line( activate_gtf.out.gtf, 30 )
    log.info "test-" + gtf_check

}


/*------------------------------------------------------------------------------------*/
/* Groovy functions
--------------------------------------------------------------------------------------*/

def boolean check_gtf_by_line( gtf_file, int line_count ) {
    def compatible = gtf_file.any { it.contains('ensembl') || it.contains('GENCODE')}
    return compatible
}

// def boolean check_gtf_by_line( f, int n ) {
//     int count = 0
//     boolean gene = false
//     boolean ensembl = false
//     boolean gencode = false

//     f.withReader('UTF-8') { r ->
//         while( count<n && ( !gene && ( !ensembl || !gencode ) ) ) {
//             line = r.readLine();
//             count = count + 1;
//             if (!gene) {
//                 if (line =~ /\bgene\b/) {
//                     gene = true
//                 }
//             };
//             if (gene && !ensembl) {
//                 if(line.contains('ensembl')) {
//                     ensembl = true
//                 }
//             };
//             if (!gene && !gencode) {
//                 if(line.contains('GENCODE')){
//                     gencode = true
//                 }
//             };
//             if (gene && ( ensembl || gencode )) {
//                 compatible = true
//             };
//         }
//     }
//     compatible
// }
