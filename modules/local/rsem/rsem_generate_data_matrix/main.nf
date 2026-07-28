process RSEM_GENERATE_DATA_MATRIX {
    label 'process_single'
    container "nciccbr/ccbr_rsem_1.3.3:v1.0"

    input:
    path gene_results     // collected *.RSEM.genes.results files
    path isoform_results  // collected *.RSEM.isoforms.results files

    output:
    path "RSEM.genes.expected_counts.${suffix}.matrix",    emit: gene_matrix
    path "RSEM.isoforms.expected_counts.${suffix}.matrix", emit: isoform_matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    suffix = task.ext.suffix ?: "all_samples"
    """
    rsem-generate-data-matrix ${gene_results.join(' ')}    > RSEM.genes.expected_counts.${suffix}.matrix
    rsem-generate-data-matrix ${isoform_results.join(' ')} > RSEM.isoforms.expected_counts.${suffix}.matrix
    """

    stub:
    suffix = task.ext.suffix ?: "all_samples"
    """
    touch RSEM.genes.expected_counts.${suffix}.matrix
    touch RSEM.isoforms.expected_counts.${suffix}.matrix
    """
}
