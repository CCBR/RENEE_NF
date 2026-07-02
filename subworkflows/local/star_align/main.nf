include { STAR_ALIGN as STAR_ALIGN_PASS1 } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as STAR_ALIGN_PASS2 } from '../../../modules/nf-core/star/align'
include { STAR_SJDB_FILTER }              from '../../../modules/local/star_sjdb_filter'

workflow star_align_workflow {
    take:
        ch_reads             // channel: [ meta, [ reads ] ]
        ch_sjdb_placeholder  // channel: path

        ch_star_index  // channel: path
        ch_star_gtf    // channel: path

    main:

        // Pass 1: align without SJDB
        STAR_ALIGN_PASS1(ch_reads, ch_star_index, ch_star_gtf, false, ch_sjdb_placeholder)

        // Filter splice junctions for Pass 2
        STAR_SJDB_FILTER(
            STAR_ALIGN_PASS1.out.spl_junc_tab
                .map { meta, sj -> sj }
                .collect()
        )

        // Pass 2: align with filtered SJDB
        STAR_ALIGN_PASS2(ch_reads, ch_star_index, ch_star_gtf, true, STAR_SJDB_FILTER.out.sjdb)

    emit:
        pass1_sj              = STAR_ALIGN_PASS1.out.spl_junc_tab
        pass1_log             = STAR_ALIGN_PASS1.out.log_final
        sjdb                  = STAR_SJDB_FILTER.out.sjdb
        pass2_log             = STAR_ALIGN_PASS2.out.log_final.mix(STAR_ALIGN_PASS2.out.log_out).mix(STAR_ALIGN_PASS2.out.log_progress)
        pass2_sj              = STAR_ALIGN_PASS2.out.spl_junc_tab
        pass2_reads_per_gene  = STAR_ALIGN_PASS2.out.read_per_gene_tab
        pass2_bam             = STAR_ALIGN_PASS2.out.bam_sorted_aligned
        pass2_transcript_bam  = STAR_ALIGN_PASS2.out.bam_transcript
}
