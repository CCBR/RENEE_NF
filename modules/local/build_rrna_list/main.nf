process BUILD_RRNA_LIST {
    label 'process_single'
    container "${params.containers.build_rnaseq}"

    input:
        path(fasta)
        path(gtf)
        val(genome_name)

    output:
        path("${genome_name}.rRNA_interval_list"), emit: rrna_list

    when:
        task.ext.when == null || task.ext.when

    script:
    // Ports RENEE (classic) workflow/rules/build.smk rule `rRNA_list`.
    """
    create_rRNA_intervals.py ${fasta} ${gtf} ${genome_name} > ${genome_name}.rRNA_interval_list
    """

    stub:
    """
    touch ${genome_name}.rRNA_interval_list
    """
}
