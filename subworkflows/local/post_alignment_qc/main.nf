include { QUALIMAP_BAMQC }              from '../../../modules/nf-core/qualimap/bamqc/main'
include { SAMTOOLS_FLAGSTAT }           from '../../../modules/CCBR/samtools/flagstat/main.nf'
include { PICARD_COLLECTRNASEQMETRICS } from '../../../modules/nf-core/picard/collectrnaseqmetrics/main.nf'
include { PRESEQ_CCURVE }               from '../../../modules/nf-core/preseq/ccurve/main'
include { HANDLE_PRESEQ_ERROR }         from '../../../modules/local/preseq/helperfunctions/main'
include { PARSE_PRESEQ_LOG }            from '../../../modules/local/preseq/helperfunctions/main'

include { RSEQC_QC } from '../rseqc_qc/main'

workflow POST_ALIGNMENT_QC {
    take:
        ch_bam         // channel: [ meta, bam ]
        ch_bai         // channel: [ meta, bai ]
        ch_bed_ref     // channel: path (genes.ref.bed)
        ch_tin_ref     // channel: path (transcripts.protein_coding_only.bed12)
        ch_refflat     // channel: path
        ch_fasta       // channel: [ meta, fasta ]
        ch_rrna_list   // channel: path
        ch_genes_gtf   // channel: [ meta, gtf ]

    main:
        ch_bam_bai = ch_bam.join(ch_bai)

        // RSeQC: strandedness, read distribution, inner distance, and TIN
        RSEQC_QC(
            ch_bam_bai,
            ch_bed_ref,
            ch_tin_ref
        )

        // Estimate library complexity and emit placeholder NRF metrics on failure
        PRESEQ_CCURVE(ch_bam)

        PRESEQ_CCURVE.out.log
            .join(ch_bam, remainder: true)
            .branch { meta, preseq_log, bam_tuple ->
                failed: preseq_log == null
                    return tuple(meta, 'nopreseqlog')
                succeeded: true
                    return tuple(meta, preseq_log)
            }
            .set { preseq_logs }

        preseq_logs.failed | HANDLE_PRESEQ_ERROR
        preseq_logs.succeeded | PARSE_PRESEQ_LOG

        ch_preseq_nrf = PARSE_PRESEQ_LOG.out.nrf
            .concat(HANDLE_PRESEQ_ERROR.out.nrf)

        // Picard RNA metrics
        PICARD_COLLECTRNASEQMETRICS(
            ch_bam,
            ch_refflat,
            ch_fasta.map { meta, fasta -> fasta },
            ch_rrna_list.ifEmpty([])
        )

        // QualiMap BAM QC
        QUALIMAP_BAMQC(
            ch_bam,
            ch_genes_gtf.map { meta, gtf -> gtf }
        )

        // Alignment summary statistics
        SAMTOOLS_FLAGSTAT(ch_bam_bai)

    emit:
        infer_experiment        = RSEQC_QC.out.infer_experiment
        read_distribution       = RSEQC_QC.out.read_distribution
        inner_distance_freq     = RSEQC_QC.out.inner_distance_freq
        inner_distance_dist     = RSEQC_QC.out.inner_distance_dist
        inner_distance_rscript  = RSEQC_QC.out.inner_distance_rscript
        tin_txt                 = RSEQC_QC.out.tin_txt
        tin_xls                 = RSEQC_QC.out.tin_xls

        preseq_ccurve           = PRESEQ_CCURVE.out.c_curve
        preseq_log              = PRESEQ_CCURVE.out.log
        preseq_nrf              = ch_preseq_nrf

        picard_rna_metrics      = PICARD_COLLECTRNASEQMETRICS.out.metrics
        picard_rna_pdf          = PICARD_COLLECTRNASEQMETRICS.out.pdf
        qualimap_results        = QUALIMAP_BAMQC.out.results
        flagstat                = SAMTOOLS_FLAGSTAT.out.flagstat
        flagstat_versions       = SAMTOOLS_FLAGSTAT.out.versions
}
