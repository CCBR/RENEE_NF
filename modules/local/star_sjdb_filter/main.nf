process STAR_SJDB_FILTER {
    label 'process_low'
    container "${params.containers.base}"
    input:
    tuple val(meta), path(sj_tabs)

    output:
    tuple val(meta), path("${meta.id}.uniq.filtered.SJ.out.tab"), emit: sjdb

    script:
    def sj_input = (sj_tabs instanceof List) ? sj_tabs.join(' ') : sj_tabs
    """
    cat $sj_input | \
        sort | \
        uniq | \
        awk -F "\t" '{if (\$5>0 && \$6==1) {print}}' | \
        cut -f1-4 | \
        sort | \
        uniq | \
        grep "^chr" | \
        grep -v "^chrM" > ${meta.id}.uniq.filtered.SJ.out.tab
    """

    stub:
    """
    touch ${meta.id}.uniq.filtered.SJ.out.tab
    """
}
