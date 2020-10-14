#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { awk } from '../../tools/luslab_linux_tools/main.nf'

workflow bt2_parse {
    take: tuple_report_meta
    main:

        // Define workflow parameters
        awk_command ='-v cols="total_reads:align1:align_gt1:non_aligned:total_aligned" \' \
        BEGIN { \
            col_count=split(cols, col_arr, ":"); \
            for(i=1; i<=col_count; i++) printf col_arr[i] ((i==col_count) ? "\\n" : ","); \
        } \
        { \
            for (i=1; i<=NF; i++) { \
                if(index($i,"reads; of these:") != 0) { \
                    split($i, line_split, " "); \
                    data["total_reads"]=line_split[1]; \
                } \
                if(index($i,"aligned concordantly exactly 1 time") != 0) { \
                    split($i, line_split, " "); \
                    data["align1"]=line_split[1]; \
                } \
                if(index($i,"aligned concordantly >1 times") != 0) { \
                    split($i, line_split, " "); \
                    data["align_gt1"]=line_split[1]; \
                } \
                if(index($i,"aligned concordantly 0 times") != 0) { \
                    split($i, line_split, " "); \
                    data["non_aligned"]=line_split[1]; \
                } \
            } \
        } \
        END { \
            data["total_aligned"] = data["align1"] + data["align_gt1"] \
            for (i=1; i<=col_count; i++) printf data[col_arr[i]] ((i==col_count) ? "\\n" : ","); \
        } \' FS="\\n" RS="\\n\\n"'

        params.modules['awk'].args = awk_command

        awk(params.modules['awk'], tuple_report_meta)

 
}