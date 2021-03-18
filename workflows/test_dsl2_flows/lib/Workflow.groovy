/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {
    static boolean check_gtf_format(gtf, max_lines) {
        def count = 0
        def line  = null
        def result = false
        gtf.withReader { reader ->
            while (line = reader.readLine()) {
                count++
                if(line.contains("GENCODE")) {
                    result = true
                    break
                }
                if(count == max_lines) {
                    break
                }     
            }
        }

        return result
    }
}
