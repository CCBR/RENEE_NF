// Helper functions to handle preseq log files

process HANDLE_PRESEQ_ERROR {
    tag "$meta.id"
    label 'process_single'
    label 'qc'
    label 'preseq'

    conda "${moduleDir}/environment.yml"
    container "${params.containers.base}"

    input:
        tuple val(meta), val(log)

    output:
        tuple val(meta), path("*nrf.txt"), emit: nrf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "NA\tNA\tNA\n" > ${prefix}.preseq.nrf.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "NA\tNA\tNA\n" > ${prefix}.preseq.nrf.txt
    """
}

process PARSE_PRESEQ_LOG {
    // Calls bin/parse_preseq_log.py to get NRF statistics from the preseq log.

    tag { meta.id }
    label 'qc'
    label 'preseq'

    container "${params.containers.base}"

    input:
        tuple val(meta), path(preseq_log)

    output:
        tuple val(meta), path("*nrf.txt"), emit: nrf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_preseq_log.py ${preseq_log} > ${prefix}.preseq.nrf.txt
    """

    stub:
    """
    touch ${meta.id}.preseqlog.nrf.txt
    """
}
