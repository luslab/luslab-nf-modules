#!/usr/bin/env nextflow
/*
========================================================================================
                         Luslab NF util module
========================================================================================
 #### Homepage / Documentation
 https://github.com/luslab/nf-modules
----------------------------------------------------------------------------------------
Module decription
----------------------------------------------------------------------------------------

This module contains helper functions for displaying help on the screen and showing 
logging output

----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.preview.dsl=2

def luslabHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """  -${c_dim}-----------------------------------------------------------------${c_reset}-
    ${c_purple} __       __    __       _______. __          ___      .______ ${c_reset}
    ${c_purple}|  |     |  |  |  |     /       ||  |        /   \\     |   _  \\ ${c_reset}
    ${c_purple}|  |     |  |  |  |    |   (----`|  |       /  ^  \\    |  |_)  | ${c_reset}
    ${c_purple}|  |     |  |  |  |     \\   \\    |  |      /  /_\\  \\   |   _  < ${c_reset}
    ${c_purple}|  `----.|  `--'  | .----)   |   |  `----./  _____  \\  |  |_)  | ${c_reset}
    ${c_purple}|_______| \\______/  |_______/    |_______/__/     \\__\\ |______/ ${c_reset}

-${c_dim}---------------------------------------------------------------${c_reset}-        
    """.stripIndent()
}