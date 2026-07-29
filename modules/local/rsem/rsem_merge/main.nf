process RSEM_MERGE {
    label 'process_single'
    container "${params.containers.base}"

    input:
    path gene_results        // collected *.RSEM.genes.results files
    path isoform_results     // collected *.RSEM.isoforms.results files
    path annotate            // annotate.genes.txt annotation file
    val  input_dir           // directory containing the results (val, used for script path)

    output:
    path "RSEM.genes.expected_count.${suffix}.txt",          emit: gene_counts
    path "RSEM.genes.FPKM.${suffix}.txt",                    emit: gene_fpkm
    path "RSEM.genes.TPM.${suffix}.txt",                     emit: gene_tpm
    path "RSEM.isoforms.expected_count.${suffix}.txt",       emit: isoform_counts
    path "RSEM.isoforms.FPKM.${suffix}.txt",                 emit: isoform_fpkm
    path "RSEM.isoforms.TPM.${suffix}.txt",                  emit: isoform_tpm
    path "RSEM.genes.expected_counts.${suffix}.reformatted.tsv", emit: reformatted

    when:
    task.ext.when == null || task.ext.when

    script:
    suffix = task.ext.suffix ?: "all_samples"
    """
    # Stage all results into a single working directory for the merge script
    merge_rsem_results.py ${annotate} . .

    if [[ "${suffix}" != "all_samples" ]]; then
        mv RSEM.genes.expected_count.all_samples.txt RSEM.genes.expected_count.${suffix}.txt
        mv RSEM.genes.FPKM.all_samples.txt RSEM.genes.FPKM.${suffix}.txt
        mv RSEM.genes.TPM.all_samples.txt RSEM.genes.TPM.${suffix}.txt
        mv RSEM.isoforms.expected_count.all_samples.txt RSEM.isoforms.expected_count.${suffix}.txt
        mv RSEM.isoforms.FPKM.all_samples.txt RSEM.isoforms.FPKM.${suffix}.txt
        mv RSEM.isoforms.TPM.all_samples.txt RSEM.isoforms.TPM.${suffix}.txt
    fi

    # Produce reformatted TSV with gene symbol
    sed 's/\\t/|/1' RSEM.genes.expected_count.${suffix}.txt | \\
        sed '1 s/^gene_id|GeneName/symbol/' > RSEM.genes.expected_counts.${suffix}.reformatted.tsv
    """

    stub:
    suffix = task.ext.suffix ?: "all_samples"
    """
    touch RSEM.genes.expected_count.${suffix}.txt
    touch RSEM.genes.FPKM.${suffix}.txt
    touch RSEM.genes.TPM.${suffix}.txt
    touch RSEM.isoforms.expected_count.${suffix}.txt
    touch RSEM.isoforms.FPKM.${suffix}.txt
    touch RSEM.isoforms.TPM.${suffix}.txt
    touch RSEM.genes.expected_counts.${suffix}.reformatted.tsv
    """
}
