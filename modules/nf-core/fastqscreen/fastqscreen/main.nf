process FASTQSCREEN_FASTQSCREEN {
    tag "${meta.id}"
    label 'process_medium'

    container "nciccbr/ccbr_fastq_screen_0.13.0:v2.0"

    input:
    tuple val(meta), path(reads)
    path fastq_screen_config
    path database // this is needed to stage the database folder to the container

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.png"), emit: png, optional: true
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    tuple val("${task.process}"), val('fastqscreen'), eval('fastq_screen --version 2>&1 | sed "s/^.*FastQ Screen v//;"'), emit: versions_fastqscreen, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    // fastq_screen names outputs after the input reads, so screening the same reads
    // against a second config collides in MultiQC unless the outputs are tagged
    def suffix = task.ext.suffix ?: ""

    """
    fastq_screen --threads ${task.cpus} \\
        --conf ${fastq_screen_config} \\
        ${reads} \\
        ${args}

    if [ -n "${suffix}" ]; then
        for f in *_screen.*; do
            mv "\$f" "\${f/_screen./${suffix}_screen.}"
        done
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_screen.html
    touch ${prefix}_screen.png
    touch ${prefix}_screen.txt
    """
}
