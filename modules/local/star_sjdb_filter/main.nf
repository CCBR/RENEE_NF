process STAR_SJDB_FILTER {
    label 'process_low'

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
        grep "^chr" | \
        grep -v "^chrM" > uniq.filtered.SJ.out.tab
    """

    stub:
    """
    touch uniq.filtered.SJ.out.tab
    """
}
