include { STAR_ALIGN as STAR_ALIGN_PASS1 } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as STAR_ALIGN_PASS2 } from '../../../modules/nf-core/star/align'
include { STAR_SJDB_FILTER }              from '../../../modules/local/star_sjdb_filter'

workflow star_align_workflow {
    take:
        ch_reads             // channel: [ meta, [ reads ] ]
        ch_star_index  // channel: path
        ch_star_gtf    // channel: path

    main:

        // Initialize a placeholder for the SJDB input to Pass 1; this will be replaced with the actual SJDB in Pass 2
        ch_sjdb_placeholder = Channel.value( file( params.sjdb_placeholder_tab, checkIfExists: true ) )

        // Pass 1: align without SJDB; combine reads with placeholder to form the input tuple
        STAR_ALIGN_PASS1(
            ch_reads.combine(ch_sjdb_placeholder),
            ch_star_index,
            ch_star_gtf,
            false
        )

        if (params.star_2_pass_basic) {
            // Run STAR_SJDB_FILTER per sample; meta flows through for SLURM-safe .join()
            STAR_SJDB_FILTER(STAR_ALIGN_PASS1.out.spl_junc_tab)

            // Join filtered SJDB into reads tuple by meta - order-safe under SLURM
            STAR_ALIGN_PASS2(
                ch_reads.join(STAR_SJDB_FILTER.out.sjdb),
                ch_star_index,
                ch_star_gtf,
                true
            )
        } else {
            // Collect all SJ tabs; use a dummy meta for the single merged STAR_SJDB_FILTER call
            STAR_SJDB_FILTER(
                STAR_ALIGN_PASS1.out.spl_junc_tab
                    .map { meta, sj -> sj }
                    .collect()
                    .map { sj_tabs -> [[id: 'merged'], sj_tabs] }
            )

            // Combine reads with the single merged SJDB (broadcast to all samples) for Pass 2
            STAR_ALIGN_PASS2(
                ch_reads.combine(STAR_SJDB_FILTER.out.sjdb.map { meta, sj -> sj }),
                ch_star_index,
                ch_star_gtf,
                true
            )
        }

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
