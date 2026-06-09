process FASTQVALIDATOR {
    tag { meta.id }
    label 'process_low'

    container "nciccbr/ccbr_fastqvalidator:v0.1.0"

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path(fastq), path('*.fastQValidator.fastq.log'), path('*.fastQValidator.exitcode.txt'), emit: result

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${fastq.baseName}"
    def minReadLen = task.ext.min_read_len ?: 2
    def logFile = "${prefix}.fastQValidator.fastq.log"
    def exitFile = "${prefix}.fastQValidator.exitcode.txt"

    """
    set +e
    fastQValidator --noeof --minReadLen ${minReadLen} --file ${fastq} > ${logFile} 2>&1
    status=\$?
    echo \$status > ${exitFile}
    exit 0
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastQValidator.fastq.log
    echo 0 > ${prefix}.fastQValidator.exitcode.txt
    """
}
