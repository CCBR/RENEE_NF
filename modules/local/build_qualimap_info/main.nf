process BUILD_QUALIMAP_INFO {
    label 'process_single'
    container "${params.containers.build_gtf}"

    input:
        path(fasta)
        path(gtf)

    output:
        path("qualimap_info.txt"), emit: qualimap_info

    when:
        task.ext.when == null || task.ext.when

    script:
    // Ports RENEE (classic) workflow/rules/build.smk rule `qualimapinfo`.
    """
    generate_qualimap_ref.py \\
        -g ${gtf} \\
        -f ${fasta} \\
        -o qualimap_info.txt \\
        --ignore-strange-chrom
    """

    stub:
    """
    touch qualimap_info.txt
    """
}
