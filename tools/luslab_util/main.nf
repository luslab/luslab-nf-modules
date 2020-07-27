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

This module contains helper functions for luslab nextflow pipelines

----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

def build_debug_param_summary() {
    Set paramsKeySet = params.keySet()
    def summary = [:]

    paramsKeySet.each {
        summary[it] = params.get(it)
    }

    output = summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    output += "\n-\033[2m----------------------------------------------------------------------\033[0m-"
    return output
}

def luslab_header() {
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

    return """  -${c_dim}----------------------------------------------------------------------${c_reset}-
    ${c_purple} __       __    __       _______. __          ___      .______ ${c_reset}
    ${c_purple}|  |     |  |  |  |     /       ||  |        /   \\     |   _  \\ ${c_reset}
    ${c_purple}|  |     |  |  |  |    |   (----`|  |       /  ^  \\    |  |_)  | ${c_reset}
    ${c_purple}|  |     |  |  |  |     \\   \\    |  |      /  /_\\  \\   |   _  < ${c_reset}
    ${c_purple}|  `----.|  `--'  | .----)   |   |  `----./  _____  \\  |  |_)  | ${c_reset}
    ${c_purple}|_______| \\______/  |_______/    |_______/__/     \\__\\ |______/ ${c_reset}

-${c_dim}----------------------------------------------------------------------${c_reset}-        

${c_green}${workflow.manifest.name} v${workflow.manifest.version}${c_reset}
    
${c_cyan}Author : ${workflow.manifest.author}${c_reset}
${c_cyan}Homepage : ${workflow.manifest.homePage}${c_reset}

-${c_dim}----------------------------------------------------------------------${c_reset}-        
    """.stripIndent()
}

// Function to ensure that resource requirements don't go beyond
// a maximum limit
def check_max(obj, type) {
  if (type == 'memory') {
    try {
      if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
        return params.max_memory as nextflow.util.MemoryUnit
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'time') {
    try {
      if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
        return params.max_time as nextflow.util.Duration
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'cpus') {
    try {
      return Math.min( obj, params.max_cpus as int )
    } catch (all) {
      println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
      return obj
    }
  }
}

def check_params(paramList) {
    Set paramsKeySet = params.keySet()

    paramList.each {
      if(!paramsKeySet.contains(it)) {
          exit 1, "Parameter " + it + " is required."
      }
      else if(params.get(it) == '') {
          exit 1, "Parameter " + it + " is required."
      }
    }
}