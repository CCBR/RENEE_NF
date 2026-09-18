process BUILD_TIN_REF {
    label 'process_single'
    container "${params.containers.build_rnaseq}"

    input:
        path(gtf)

    output:
        path("transcripts.protein_coding_only.bed12"), emit: tin_ref

    when:
        task.ext.when == null || task.ext.when

    script:
    // Ports RENEE (classic) workflow/rules/build.smk rule `tin_ref`.
    """
    gtf2protein_coding_genes.py ${gtf} > protein_coding_genes.lst
    gtfToGenePred -ignoreGroupsWithoutExons -genePredExt ${gtf} genes.gtf.genePred
    genePredToBed genes.gtf.genePred genes.gtf.genePred.bed

    awk -F '\\t' -v OFS='\\t' '{print \$12,\$1}' genes.gtf.genePred \\
        | sort -k1,1n > gene2transcripts

    while read gene; do
        awk -F '\\t' -v gene="\${gene}" '\$1 == gene' gene2transcripts;
    done < protein_coding_genes.lst > gene2transcripts.protein_coding_only

    gene2transcripts_add_length.py \\
        gene2transcripts.protein_coding_only \\
        genes.gtf.genePred.bed \\
        > gene2transcripts.protein_coding_only.with_len

    sort -k1,1 -k3,3nr gene2transcripts.protein_coding_only.with_len | \\
        awk -F '\\t' '{if (!seen[\$1]) {seen[\$1]++; print \$2}}' > protein_coding_only.txt

    # protein_coding_only.txt's transcript IDs are stripped of their version
    # suffix (the awk -F '.' below) before matching, since versions can churn
    # between GTF releases -- so field 4 (name) is stripped the same way
    # here for an exact, version-insensitive comparison, rather than relying
    # on the stripped ID being an (unanchored, collision-prone) substring of
    # the full versioned line.
    while read transcript; do
        awk -F '\\t' -v transcript="\${transcript}" '{split(\$4, base, "."); if (base[1] == transcript) {print; exit}}' genes.gtf.genePred.bed;
    done < <(awk -F '.' '{print \$1}' protein_coding_only.txt) > transcripts.protein_coding_only.bed12

    rm -f "protein_coding_genes.lst" "genes.gtf.genePred" "genes.gtf.genePred.bed" \\
        "gene2transcripts" "gene2transcripts.protein_coding_only" "protein_coding_only.txt" \\
        "gene2transcripts.protein_coding_only.with_len"
    """

    stub:
    """
    touch transcripts.protein_coding_only.bed12
    """
}
