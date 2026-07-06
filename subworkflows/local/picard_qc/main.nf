// CRAFT: Port Snakemake Picard post-alignment QC while reusing shared modules.

include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main.nf'
include { PICARD_MARKDUPLICATES }         from '../../../modules/nf-core/picard/markduplicates/main.nf'
include { PICARD_COLLECTRNASEQMETRICS }   from '../../../modules/nf-core/picard/collectrnaseqmetrics/main.nf'
include { SAMTOOLS_FLAGSTAT }             from '../../../modules/CCBR/samtools/flagstat/main.nf'

process RENEE_NORMALIZE_PICARD_DUPLICATE_METRICS {
    tag "${meta.id}"
    label 'process_single'

    container "${params.containers.base}"

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path("*.star.duplic"), emit: metrics

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/MarkDuplicates/picard.sam.MarkDuplicates/g' \\
        ${metrics} \\
        > ${prefix}.star.duplic
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.star.duplic
    """
}

process RENEE_NORMALIZE_PICARD_RNASEQ_METRICS {
    tag "${meta.id}"
    label 'process_single'

    container "${params.containers.base}"

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path("*.RnaSeqMetrics.txt"), emit: metrics

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' \\
        ${metrics} \\
        > ${prefix}.RnaSeqMetrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.RnaSeqMetrics.txt
    """
}

process RENEE_FLAGSTAT_CONCORD {
    tag "${meta.id}"
    label 'process_single'

    container "${params.containers.base}"

    input:
    tuple val(meta), path(flagstat)

    output:
    tuple val(meta), path("*.flagstat.concord.txt"), emit: flagstat_concord

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${flagstat} > ${prefix}.flagstat.concord.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.flagstat.concord.txt
    """
}

workflow PICARD_QC {
    take:
        ch_bam       // channel: [ meta, bam ]
        ch_fasta     // channel: path
        ch_refflat   // channel: path
        ch_rrna_list // channel: path

    main:
        ch_empty_reference = Channel.value([ [:], [], [] ])

        PICARD_ADDORREPLACEREADGROUPS(ch_bam, ch_empty_reference)
        PICARD_MARKDUPLICATES(PICARD_ADDORREPLACEREADGROUPS.out.bam, ch_empty_reference)

        ch_markdup_bam = PICARD_MARKDUPLICATES.out.bam
        ch_markdup_bai = PICARD_MARKDUPLICATES.out.bai

        PICARD_COLLECTRNASEQMETRICS(ch_markdup_bam, ch_refflat, ch_fasta, ch_rrna_list)

        SAMTOOLS_FLAGSTAT(ch_markdup_bam.join(ch_markdup_bai))

        RENEE_NORMALIZE_PICARD_DUPLICATE_METRICS(PICARD_MARKDUPLICATES.out.metrics)
        RENEE_NORMALIZE_PICARD_RNASEQ_METRICS(PICARD_COLLECTRNASEQMETRICS.out.metrics)
        RENEE_FLAGSTAT_CONCORD(SAMTOOLS_FLAGSTAT.out.flagstat)

    emit:
        bam               = ch_markdup_bam
        bai               = ch_markdup_bai
        duplicate_metrics = RENEE_NORMALIZE_PICARD_DUPLICATE_METRICS.out.metrics
        rnaseq_metrics    = RENEE_NORMALIZE_PICARD_RNASEQ_METRICS.out.metrics
        flagstat_concord  = RENEE_FLAGSTAT_CONCORD.out.flagstat_concord
}
