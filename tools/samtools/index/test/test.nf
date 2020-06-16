#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting Samtools index module")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include samtools_index from '../main.nf' //addParams(star_custom_args: "-t 2")

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testMetaData = [
  ['Sample1', "$baseDir/input/prpf8_eif4a3_rep1.Unmapped.fqAligned.sortedByCoord.out.bam"],
  ['Sample2', "$baseDir/input/prpf8_eif4a3_rep2.Unmapped.fqAligned.sortedByCoord.out.bam"]
]

 Channel
    .from( testMetaData )
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set { ch_testData }

/*------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // Run star
    samtools_index( ch_testData )

    // Collect file names and view output
    samtools_index.out.baiFiles | view
}
