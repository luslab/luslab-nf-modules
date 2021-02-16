#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Processes
--------------------------------------------------------------------------------------*/

process count_lines {
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    tag "${sample_id}"

    container 'ubuntu:16.04'

    input:
        tuple val(sample_id), path(input_file)

    output:
        tuple val(sample_id), stdout, emit: line_count

    // Treats gzipped and uncompressed files similarly
    script:
    """
    echo -n "\$(zcat -f $input_file | wc -l)"
    """
}

/*------------------------------------------------------------------------------------*/
/* Workflows
--------------------------------------------------------------------------------------*/

// Check number of items in an output channel
workflow assert_channel_count {
    take: test_channel
    take: channel_name
    take: expected
    main:
        test_channel.count()
            .subscribe{
                channel_count = it
                if(channel_count != expected) {
                    throw new Exception(channel_name + " channel count is " + channel_count + ", expected count is " + expected);
                }
        }
}

workflow assert_line_number {
    take: test_channel
    take: channel_name
    take: expected_line_counts
    main:
        count_lines(test_channel)

        count_lines.out.subscribe {
            if(expected_line_counts[it[0]] != it[1].toInteger()) {
                throw new Exception("Error with channel " + channel_name + ": Sample " + it[0] + " is expected to have " + expected_line_counts[it[0]] + " lines, but has " + it[1])
            }
        }
}
