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

process mn_cpu_resources {
    label 'mn_cpu'
    script:
      message = "mn_cpu - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process l_cpu_resources {
    label 'l_cpu'
    script:
      message = "l_cpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process m_cpu_resources {
    label 'm_cpu'
    script:
      message = "m_cpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process h_cpu_resources {
    label 'h_cpu'
    script:
      message = "h_cpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process mx_cpu_resources {
    label 'mx_cpu'
    script:
      message = "mx_cpu - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process mn_hmem_resources {
    label 'mn_hmem'
    script:
      message = "mn_hmem - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process l_hmem_resources {
    label 'l_hmem'
    script:
      message = "l_hmem  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process m_hmem_resources {
    label 'm_hmem'
    script:
      message = "m_hmem  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process h_hmem_resources {
    label 'h_hmem'
    script:
      message = "h_hmem  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process mx_hmem_resources {
    label 'mx_hmem'
    script:
      message = "mx_hmem - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process mn_gpu_resources {
    label 'mn_gpu'
    script:
      message = "mn_gpu - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process l_gpu_resources {
    label 'l_gpu'
    script:
      message = "l_gpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process m_gpu_resources {
    label 'm_gpu'
    script:
      message = "m_gpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process h_gpu_resources {
    label 'h_gpu'
    script:
      message = "h_gpu  - cpus=" + task.cpus
      message += " mem=" + task.memory.toString().replace(' ', '')
      log.info message
      """
      echo ${message}
      """
}

process mx_gpu_resources {
    label 'mx_gpu'
    script:
      message = "mx_gpu - cpus=" + task.cpus
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
        mn_cpu_resources()
        l_cpu_resources()
        m_cpu_resources()
        h_cpu_resources()
        mx_cpu_resources()
        mn_hmem_resources()
        l_hmem_resources()
        m_hmem_resources()
        h_hmem_resources()
        mx_hmem_resources()
        mn_gpu_resources()
        l_gpu_resources()
        m_gpu_resources()
        h_gpu_resources()
        mx_gpu_resources()

        // Doesnt really work as we can only access the cpu and mem variables. Have to resort to manual testing on the cluster
        // default_expected_cpu = 2
        // default_resources.out.cpus.subscribe{
        //         if(it != default_expected_cpu) {
        //             throw new Exception("Default resources CPU count is " + it + " when is should be " + default_expected_cpu);
        //         }
        // }
}