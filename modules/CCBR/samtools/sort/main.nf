process SAMTOOLS_SORT {
    tag { meta.id }
    label 'process_high'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    tuple val(meta), path("*.csi"),                emit: csi, optional: true
    path("versions.yml"),                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bam.baseName}.sorted"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # current working directory is a tmpdir when 'scratch' is set
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT

    samtools sort \\
        ${args} \\
        -@ ${task.cpus} \\
        -o ${prefix}.bam##idx##${prefix}.bam.bai \\
        -m 2G \\
        -T \$TMP \\
        --write-index \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${bam.baseName}.sorted"
    """
    touch ${prefix}.bam ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
