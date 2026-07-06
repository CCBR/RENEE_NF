include { STAR_ALIGN as STAR_ALIGN_PASS1  } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as STAR_ALIGN_PASS2  } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as STAR_ALIGN_BASIC  } from '../../../modules/nf-core/star/align'
include { STAR_SJDB_FILTER }               from '../../../modules/local/star_sjdb_filter'

workflow STAR_ALIGN {
    take:
        ch_reads       // channel: [ meta, [ reads ] ]
        ch_star_index  // channel: path
        ch_star_gtf    // channel: path

    main:
        ch_sjdb_placeholder = Channel.value(file(params.sjdb_placeholder_tab, checkIfExists: true))

        if (params.star_2_pass_basic) {
            // Per-sample two-pass: STAR handles both passes internally via --twopassMode Basic
            STAR_ALIGN_BASIC(ch_reads, ch_star_index, ch_star_gtf, false, ch_sjdb_placeholder)

            ch_pass1_sj             = Channel.empty()
            ch_pass1_log            = Channel.empty()
            ch_sjdb                 = Channel.empty()
            ch_pass2_log            = STAR_ALIGN_BASIC.out.log_final
                                          .mix(STAR_ALIGN_BASIC.out.log_out)
                                          .mix(STAR_ALIGN_BASIC.out.log_progress)
            ch_pass2_sj             = STAR_ALIGN_BASIC.out.spl_junc_tab
            ch_pass2_reads_per_gene = STAR_ALIGN_BASIC.out.read_per_gene_tab
            ch_pass2_bam            = STAR_ALIGN_BASIC.out.bam_sorted_aligned
            ch_pass2_transcript_bam = STAR_ALIGN_BASIC.out.bam_transcript
        } else {
            // Multi-sample two-pass alignment
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

            ch_pass1_sj             = STAR_ALIGN_PASS1.out.spl_junc_tab
            ch_pass1_log            = STAR_ALIGN_PASS1.out.log_final
            ch_sjdb                 = STAR_SJDB_FILTER.out.sjdb
            ch_pass2_log            = STAR_ALIGN_PASS2.out.log_final
                                          .mix(STAR_ALIGN_PASS2.out.log_out)
                                          .mix(STAR_ALIGN_PASS2.out.log_progress)
            ch_pass2_sj             = STAR_ALIGN_PASS2.out.spl_junc_tab
            ch_pass2_reads_per_gene = STAR_ALIGN_PASS2.out.read_per_gene_tab
            ch_pass2_bam            = STAR_ALIGN_PASS2.out.bam_sorted_aligned
            ch_pass2_transcript_bam = STAR_ALIGN_PASS2.out.bam_transcript
        }

    emit:
        pass1_sj              = ch_pass1_sj
        pass1_log             = ch_pass1_log
        sjdb                  = ch_sjdb
        pass2_log             = ch_pass2_log
        pass2_sj              = ch_pass2_sj
        pass2_reads_per_gene  = ch_pass2_reads_per_gene
        pass2_bam             = ch_pass2_bam
        pass2_transcript_bam  = ch_pass2_transcript_bam
}
