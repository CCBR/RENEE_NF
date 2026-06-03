process CUTADAPT {
    tag { meta.id }
    label 'process_high'

    container 'nciccbr/ncigb_cutadapt_v1.18:latest'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*.trim.fastq.gz'), emit: reads
        tuple val(meta), path('*.log')          , emit: log
        path "versions.yml"                     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    def args = [
            '--nextseq-trim=2',
            '--trim-n -n 5 -O 5',
            '-q 10,10',
            '-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa'
        ]
    if (meta.single_end) {
        args += [
            '-m 20'
        ]
    } else {
        args += [
            '-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa',
            '-m 20:20',
        ]
    }
    args = args.join(' ').trim()
    """
    cutadapt \\
        --cores ${task.cpus} \\
        ${args} \\
        ${trimmed} \\
        ${reads} \\
        > ${prefix}.cutadapt.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trim.fastq.gz" : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
