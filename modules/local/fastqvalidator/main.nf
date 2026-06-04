process FASTQVALIDATOR {
    tag { meta.id }
    label 'process_low'

    container "${params.containers.fastqvalidator}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*.fastQValidator.txt'), emit: report
        tuple val(meta), path('*.fastQValidator.R1.fastq.log'), emit: log_r1
        tuple val(meta), path('*.fastQValidator.R2.fastq.log'), optional: true, emit: log_r2

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minReadLen = task.ext.min_read_len ?: 2
    def readList = reads instanceof List ? reads : [reads]
    def r1 = readList[0]
    def r2 = readList.size() > 1 ? readList[1] : null
    def reportFile = "${prefix}.fastQValidator.txt"
    def r1Log = "${prefix}.fastQValidator.R1.fastq.log"
    def r2Log = "${prefix}.fastQValidator.R2.fastq.log"
    def pairedCommands = r2 ? """
    fastQValidator --noeof --minReadLen ${minReadLen} --file ${r2} > ${r2Log}
    cat ${r1Log} ${r2Log} > ${reportFile}
    """ : """
    cp ${r1Log} ${reportFile}
    """
    """
    fastQValidator --noeof --minReadLen ${minReadLen} --file ${r1} > ${r1Log}
    ${pairedCommands}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastQValidator.R1.fastq.log
    touch ${prefix}.fastQValidator.R2.fastq.log
    touch ${prefix}.fastQValidator.txt
    """
}
