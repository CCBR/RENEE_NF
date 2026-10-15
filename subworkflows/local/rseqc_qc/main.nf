include { RSEQC_INFEREXPERIMENT  } from '../../../modules/nf-core/rseqc/inferexperiment/main.nf'
include { RSEQC_READDISTRIBUTION } from '../../../modules/nf-core/rseqc/readdistribution/main.nf'
include { RSEQC_INNERDISTANCE    } from '../../../modules/nf-core/rseqc/innerdistance/main.nf'
include { RSEQC_TIN              } from '../../../modules/nf-core/rseqc/tin/main.nf'

workflow RSEQC_QC {
    take:
        ch_bam_bai   // channel: [ meta, bam, bai ]
        ch_bed_ref   // channel: path  (genes.ref.bed)
        ch_tin_ref   // channel: path  (transcripts.protein_coding_only.bed12)

    main:
        // Infer strandedness from read orientation
        RSEQC_INFEREXPERIMENT(ch_bam_bai, ch_bed_ref)

        // Compute read distributions over genomic features
        RSEQC_READDISTRIBUTION(ch_bam_bai, ch_bed_ref)

        // Calculate inner distance between paired-end read mates.
        // inner_distance.py requires mate pairs, so single-end samples are
        // filtered out (RENEE's Snakemake inner_distance rule is paired-end only).
        RSEQC_INNERDISTANCE(
            ch_bam_bai.filter { meta, bam, bai -> !meta.single_end },
            ch_bed_ref
        )

        // Compute transcript integrity numbers (TIN) for canonical protein-coding genes
        RSEQC_TIN(ch_bam_bai, ch_tin_ref)

    emit:
        infer_experiment        = RSEQC_INFEREXPERIMENT.out.txt        // [ meta, *.infer_experiment.txt ]
        read_distribution       = RSEQC_READDISTRIBUTION.out.txt        // [ meta, *.read_distribution.txt ]
        inner_distance_freq     = RSEQC_INNERDISTANCE.out.freq          // [ meta, *freq.txt ]
        inner_distance_dist     = RSEQC_INNERDISTANCE.out.distance      // [ meta, *distance.txt ]
        inner_distance_rscript  = RSEQC_INNERDISTANCE.out.rscript       // [ meta, *.r ]
        tin_txt                 = RSEQC_TIN.out.txt                     // [ meta, *.txt ]
        tin_xls                 = RSEQC_TIN.out.xls                     // [ meta, *.xls ]
}
