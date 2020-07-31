#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Check number of items in an output channel
workflow assert_channel_count {
    take: test_channel
    take: channelName
    take: expected
    main:
        test_channel.count()
            .subscribe{
                channelCount = it
                if(channelCount != expected) {
                    throw new Exception(channelName + " channel count is " + channelCount + ", expected count is " + expected);
                }
        }
}