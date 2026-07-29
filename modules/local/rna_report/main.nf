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

    TIN_FILES=( ${tin_files.join(' ')} )
    if [[ \${#TIN_FILES[@]} -eq 1 ]] && [[ "\$(head -n 1 "\${TIN_FILES[0]}")" != \$'geneID\\tchrom\\ttx_start\\ttx_end\\tTIN' ]]; then
        cp "\${TIN_FILES[0]}" combined_TIN.tsv
    else
        Rscript - "\${TIN_FILES[@]}" <<'RSCRIPT'
    inputs <- commandArgs(trailingOnly = TRUE)
    tin_tables <- lapply(inputs, function(input) {
        tin <- read.delim(input, check.names = FALSE)
        sample_name <- sub("\\\\.tin\\\\.xls\$", "", basename(input))
        sample_tin <- data.frame(
            transcript_id = tin[["geneID"]],
            tin = tin[["TIN"]],
            check.names = FALSE
        )
        names(sample_tin)[2] <- sample_name
        sample_tin
    })
    combined_tin <- Reduce(
        function(left, right) merge(left, right, by = "transcript_id", all = TRUE),
        tin_tables
    )
    write.table(
        combined_tin,
        "combined_TIN.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )
    RSCRIPT
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
