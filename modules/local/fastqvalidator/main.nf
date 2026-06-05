process FASTQVALIDATOR {
    tag { meta.id }
    label 'process_low'

    container "nciccbr/ccbr_fastqvalidator:v0.1.0"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*.fastQValidator.txt'), emit: report
        tuple val(meta), path('*.fastQValidator.R1.fastq.log'), emit: log_r1
        tuple val(meta), path('*.fastQValidator.R2.fastq.log'), optional: true, emit: log_r2
        tuple val(meta), path('*.fastQValidator.status.txt'), emit: status

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
    def statusFile = "${prefix}.fastQValidator.status.txt"

    // Run both validations and capture their statuses so workflow logic can decide whether to stop.
    def pairedCommands = r2 ? """
    r2_status=0
    fastQValidator --noeof --minReadLen ${minReadLen} --file ${r2} > ${r2Log} 2>&1 || r2_status=\$?
    cat ${r1Log} ${r2Log} > ${reportFile}
    """ : """
    r2_status=0
    cp ${r1Log} ${reportFile}
    """
    """
    r1_status=0
    fastQValidator --noeof --minReadLen ${minReadLen} --file ${r1} > ${r1Log} 2>&1 || r1_status=\$?
    ${pairedCommands}

    if [[ "\${r1_status}" -eq 0 && "\${r2_status}" -eq 0 ]]; then
        echo "PASS" > ${statusFile}
    else
        echo "FAIL\tr1=\${r1_status}\tr2=\${r2_status}" > ${statusFile}
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastQValidator.R1.fastq.log
    touch ${prefix}.fastQValidator.R2.fastq.log
    touch ${prefix}.fastQValidator.txt
    echo "PASS" > ${prefix}.fastQValidator.status.txt
    """
}
