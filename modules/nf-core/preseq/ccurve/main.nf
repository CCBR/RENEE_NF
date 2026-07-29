process PRESEQ_CCURVE {
    tag "$meta.id"
    label 'process_single'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "nciccbr/ccbr_preseq:v0.0.1"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.c_curve.txt"), emit: c_curve
    tuple val(meta), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    preseq \\
        c_curve \\
        $args \\
        -output ${prefix}.c_curve.txt \\
        $bam \\
        2> ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.c_curve.txt
    touch ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
