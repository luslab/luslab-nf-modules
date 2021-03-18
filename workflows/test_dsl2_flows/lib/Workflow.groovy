/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {
    static void check_gtf_format(gtf) {
        def count = 0
        def line  = null
        gtf.withReader { reader ->
            while (line = reader.readLine()) {
            }
        }
    }
}
