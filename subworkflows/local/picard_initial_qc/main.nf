include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main.nf'
include { PICARD_MARKDUPLICATES }         from '../../../modules/nf-core/picard/markduplicates/main.nf'

workflow PICARD_INITIAL_QC {
    take:
        ch_bam       // channel: [ meta, bam ]
    main:

        // Renee does not use reference fastas, this is a dummy ch for it
        ch_empty_reference = Channel.value([ [:], [], [] ])

        PICARD_ADDORREPLACEREADGROUPS(ch_bam, ch_empty_reference)
        PICARD_MARKDUPLICATES(PICARD_ADDORREPLACEREADGROUPS.out.bam, ch_empty_reference)

        ch_markdup_bam = PICARD_MARKDUPLICATES.out.bam
        ch_markdup_bai = PICARD_MARKDUPLICATES.out.bai

    emit:
        bam = PICARD_MARKDUPLICATES.out.bam // [ meta, *.bam ]
        bai = PICARD_MARKDUPLICATES.out.bai // [ meta, *.bai ]
}
