#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

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