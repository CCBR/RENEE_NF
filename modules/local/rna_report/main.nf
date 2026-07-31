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
    def tin_files = tins instanceof List ? tins : [tins]
    """
    export R_LIBS_SITE=/usr/local/lib/R/site-library

    # Skip TIN matrix creation if only one TIN file is provided and it has the correct header

    TIN_FILES=( ${tin_files.join(' ')} )
    if [[ \${#TIN_FILES[@]} -eq 1 ]] && [[ "\$(head -n 1 "\${TIN_FILES[0]}")" != \$'geneID\\tchrom\\ttx_start\\ttx_end\\tTIN' ]]; then
        cp "\${TIN_FILES[0]}" combined_TIN.tsv
    else
        create_tin_matrix.py "\${TIN_FILES[@]}" > combined_TIN.tsv
    fi

    rna_report.R \
        -m ${rmarkdown} \
        -r "\$PWD/${counts}" \
        -t "\$PWD/combined_TIN.tsv" \
        -q "\$PWD/${qc}" \
        -f ${prefix}.html \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'RNA_Report'
    """
    touch ${prefix}.html
    """
}
