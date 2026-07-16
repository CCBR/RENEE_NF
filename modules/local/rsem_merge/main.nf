process RSEM_MERGE {
    label 'process_single'
    container "${params.containers.base}"

    input:
    path gene_results        // collected *.RSEM.genes.results files
    path isoform_results     // collected *.RSEM.isoforms.results files
    path annotate            // annotate.genes.txt annotation file
    val  input_dir           // directory containing the results (val, used for script path)

    output:
    path "RSEM.genes.expected_count.all_samples.txt",          emit: gene_counts
    path "RSEM.genes.FPKM.all_samples.txt",                    emit: gene_fpkm
    path "RSEM.genes.TPM.all_samples.txt",                     emit: gene_tpm
    path "RSEM.isoforms.expected_count.all_samples.txt",       emit: isoform_counts
    path "RSEM.isoforms.FPKM.all_samples.txt",                 emit: isoform_fpkm
    path "RSEM.isoforms.TPM.all_samples.txt",                  emit: isoform_tpm
    path "RSEM.genes.expected_counts.all_samples.reformatted.tsv", emit: reformatted

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Stage all results into a single working directory for the merge script
    merge_rsem_results.py ${annotate} . .

    # Produce reformatted TSV with gene symbol
    sed 's/\\t/|/1' RSEM.genes.expected_count.all_samples.txt | \\
        sed '1 s/^gene_id|GeneName/symbol/' > RSEM.genes.expected_counts.all_samples.reformatted.tsv
    """

    stub:
    """
    touch RSEM.genes.expected_count.all_samples.txt
    touch RSEM.genes.FPKM.all_samples.txt
    touch RSEM.genes.TPM.all_samples.txt
    touch RSEM.isoforms.expected_count.all_samples.txt
    touch RSEM.isoforms.FPKM.all_samples.txt
    touch RSEM.isoforms.TPM.all_samples.txt
    touch RSEM.genes.expected_counts.all_samples.reformatted.tsv
    """
}
