#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Processes
--------------------------------------------------------------------------------------*/

process default_resources {
    script:
      message = "default   - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process max_cpu_q_resources {
    script:
      message = "max_cpu_q - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process cpu_q_mncpu_mnmem_resources {
    label 'mncpu_mnmem'
    script:
      message = "cpu_q_mn-cpu_mn-mem - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

/*------------------------------------------------------------------------------------*/
/* Workflows
--------------------------------------------------------------------------------------*/

// Check the resource allocation models for different profiles and configurations
workflow assert_resource_allocation_models {
    main:
        default_resources()
        max_cpu_q_resources()
        cpu_q_mncpu_mnmem_resources()

        // Doesnt really work as we can only access the cpu and mem variables. Have to resort to manual testing on the cluster
        // default_expected_cpu = 2
        // default_resources.out.cpus.subscribe{
        //         if(it != default_expected_cpu) {
        //             throw new Exception("Default resources CPU count is " + it + " when is should be " + default_expected_cpu);
        //         }
        // }
}