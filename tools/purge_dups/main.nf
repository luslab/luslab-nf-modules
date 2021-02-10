#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process purge_dups {
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'quay.io/biocontainers/purge_dups:1.2.5--hed695b0_0'

    input:
        val opts
        tuple val(meta), path(primary_assembly)
        tuple val(meta), path(aligned_reads)

    output:
        tuple val(meta), path("purged.fa"), emit: purged_fasta
        tuple val(meta), path("hap.fa"), emit: hap_fasta

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        // See https://github.com/dfguan/purge_dups#--pipeline-guide
        pbcstat_command    = "/usr/local/bin/pbcstat ${aligned_reads}"  // produces PB.base.cov and PB.stat files
        calcuts_command    = "/usr/local/bin/calcuts PB.stat > cutoffs 2> calcults.log"
        split_fa_command   = "/usr/local/bin/split_fa ${primary_assembly} > ${primary_assembly}.split"
        minimap_command    = "/usr/local/bin/minimap2 -xasm5 -DP ${primary_assembly}.split ${primary_assembly}.split | gzip > ${primary_assembly}.split.self.paf.gz"
        purge_dups_command = "/usr/local/bin/purge_dups -2 -T cutoffs -c PB.base.cov ${primary_assembly}.split.self.paf.gz > dups.bed 2> purge_dups.log"
        get_seqs_commmand  = "/usr/local/bin/get_seqs -e dups.bed ${primary_assembly}"

        if (params.verbose){
            println ("[MODULE] purge_dups calculate coverage command: " + pbcstat_command)
            println ("[MODULE] purge_dups calculate cutoff command: " + calcuts_command)
            println ("[MODULE] purge_dups calculate split FASTA command: " + split_fa_command)
            println ("[MODULE] purge_dups calculate self-align command: " + minimap_command)
            println ("[MODULE] purge_dups calculate purge dups command: " + purge_dups_command)
            println ("[MODULE] purge_dups calculate get sequences command: " + get_seqs_commmand)
        }

    """
    ${pbcstat_command}
    ${calcuts_command}
    ${split_fa_command}
    ${minimap_command}
    ${purge_dups_command}
    ${get_seqs_commmand}
    """
}
