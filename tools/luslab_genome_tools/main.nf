#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Converts region in the format - "chrX:15560138-15602945 to a bed file for downstream use"
process region2bed {
    input:
        tuple val(region)

    output:
        path "custom_region.bed", emit: bedFile

    script:

    // Log
    if (params.verbose){
        println ("[MODULE] region2bed region: " + region)
    }
    
    //SHELL
    """
    echo $region | awk '{ \
    split(\$0,a,":"); \
    split(a[2],b,"-"); \
    printf "%s\t%s\t%s",a[1],b[1],b[2]}' \
    > custom_region.bed
    """
}

