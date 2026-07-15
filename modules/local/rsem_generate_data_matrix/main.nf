process RSEM_GENERATE_DATA_MATRIX {
    label 'process_single'
    container "nciccbr/ccbr_rsem_1.3.3:v1.0"

    input:
    path gene_results     // collected *.RSEM.genes.results files
    path isoform_results  // collected *.RSEM.isoforms.results files

    output:
    path "RSEM.genes.expected_counts.all_samples.matrix",    emit: gene_matrix
    path "RSEM.isoforms.expected_counts.all_samples.matrix", emit: isoform_matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rsem-generate-data-matrix ${gene_results.join(' ')}    > RSEM.genes.expected_counts.all_samples.matrix
    rsem-generate-data-matrix ${isoform_results.join(' ')} > RSEM.isoforms.expected_counts.all_samples.matrix
    """

    stub:
    """
    touch RSEM.genes.expected_counts.all_samples.matrix
    touch RSEM.isoforms.expected_counts.all_samples.matrix
    """
}
