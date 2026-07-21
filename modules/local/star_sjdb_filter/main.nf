process STAR_SJDB_FILTER {
    label 'process_low'
    container "${params.containers.base}"
    input:
    path sj_tabs

    output:
    path 'uniq.filtered.SJ.out.tab', emit: sjdb

    script:
    """
    cat ${sj_tabs.join(' ')} | \
        sort | \
        uniq | \
        awk -F "\t" '{if (\$5>0 && \$6==1) {print}}' | \
        cut -f1-4 | \
        sort | \
        uniq | \
        awk '/^chr/ && !/^chrM/' > uniq.filtered.SJ.out.tab
        # awk is used instead of grep to avoid a non-zero exit code when no
        # chr-prefixed junctions are present (e.g. sarscov2 test data). grep
        # exits 1 on no match; awk always exits 0, so pipeline errors are real.
    """

    stub:
    """
    touch uniq.filtered.SJ.out.tab
    """
}
