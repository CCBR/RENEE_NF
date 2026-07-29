process RNA_REPORT {
    tag "$meta.id"
    label 'process_low'

    container "${params.containers.rna}"

    input:
    tuple val(meta), path(counts), path(tins), path(qc), path(rmarkdown)

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val("${task.process}"), val('R'), eval("R --version | head -n 1 | cut -d ' ' -f 3"), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'RNA_Report'
    """
    export R_LIBS_SITE=/usr/local/lib/R/site-library

    rna_report.R \
        -m ${rmarkdown} \
        -r ${counts} \
        -t ${tins} \
        -q ${qc} \
        -f ${prefix}.html \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'RNA_Report'
    """
    touch ${prefix}.html
    """
}
