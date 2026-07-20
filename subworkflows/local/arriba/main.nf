include { FUSION_ARRIBA                             } from '../../../modules/local/fusion_arriba/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ARRIBA    } from '../../../modules/CCBR/samtools/sort/main.nf'
include { ARRIBA_VISUALISATION                     } from '../../../modules/nf-core/arriba/visualisation/main.nf'

// Arriba gene-fusion calling subworkflow.
//
// Mirrors the paired-end `arriba` Snakemake rule from RENEE:
//   1. STAR alignment with Arriba chimeric-detection parameters
//      (--twopassMode Basic, all --chim* flags — see modules.config STAR_ALIGN_ARRIBA block).
//   2. Arriba fusion calling on the unsorted chimeric BAM.
//   3. samtools sort + index of the chimeric BAM (required for visualisation).
//   4. draw_fusions.R figure generation via ARRIBA_VISUALISATION.
//
// The workflow is only called when the genome supplies a blacklist file.
// cytobands and protein_domains are optional; visualisation degrades gracefully
// when they are absent (flags are omitted from the draw_fusions.R call).

workflow arriba {

    take:
        ch_reads            // channel: [ val(meta), path(reads) ]         trimmed FASTQ reads (paired-end)
        ch_star_index       // channel: [ val(meta2), path(index) ]        STAR genome index
        ch_genes_gtf        // channel: [ val(meta3), path(gtf) ]          annotation GTF
        ch_fasta            // channel: [ val(meta), path(fasta) ]          genome FASTA
        ch_blacklist        // channel: path(blacklist)     — may be empty
        ch_known_fusions    // channel: path(known_fusions) — may be empty
        ch_cytobands        // channel: path(cytobands)     — may be empty
        ch_protein_domains  // channel: path(protein_domains) — may be empty

    main:

        // Scope note: Arriba is executed via a single Snakemake rule in
        // `workflow/rules/paired-end.smk`. The rule is conditional — it only runs when the
        // three required reference files (`FUSIONBLACKLIST`, `FUSIONCYTOBAND`,
        // `FUSIONPROTDOMAIN`) are present in the genome config. The Arriba container image
        // is also shared with all STAR alignment rules. Single-end mode does not have a
        // dedicated Arriba fusion-calling rule.

        // Gate on the blacklist: if the genome config does not supply a fusion blacklist
        // (ch_blacklist is empty), combine() produces no items and all downstream
        // processes are skipped — mirroring the Snakemake conditional rule.
        // Single-end reads are also excluded: the Snakemake arriba rule is paired-end only.
        //
        // The blacklist is carried through combine() so ch_blacklist is only consumed once
        // (queue channels cannot be safely reused across multiple operators).
        ch_reads
            .filter { meta, reads -> !meta.single_end }
            .combine( ch_blacklist )
            .multiMap { meta, reads, bl ->
                reads_in:     tuple( meta, reads )
                blacklist_in: bl
            }
            .set { ch_fusion_arriba_inputs }

        FUSION_ARRIBA(
            ch_fusion_arriba_inputs.reads_in,
            ch_fasta,
            ch_genes_gtf,
            ch_star_index,
            ch_fusion_arriba_inputs.blacklist_in
        )

        // Sort and index the chimeric BAM (CCBR samtools/sort writes BAI in one step).
        // Prefix is set to "${meta.id}.arriba" via modules.config SAMTOOLS_SORT_ARRIBA block.
        SAMTOOLS_SORT_ARRIBA( FUSION_ARRIBA.out.bam )

        // Join sorted BAM + BAI with fusions TSV → [ meta, bam, bai, fusions ]
        // to satisfy the ARRIBA_VISUALISATION input tuple.
        ch_vis_input = SAMTOOLS_SORT_ARRIBA.out.bam
            .join( FUSION_ARRIBA.out.fusions )

        // Wrap optional reference paths as tuples for ARRIBA_VISUALISATION.
        // Falls back to a null-meta empty tuple so the process flag is omitted
        // when the reference file is absent.
        ch_prot_vis = ch_protein_domains
            .map { pd -> [ [id: 'reference'], pd ] }
            .ifEmpty( [ [id: 'null'], [] ] )

        ch_cyto_vis = ch_cytobands
            .map { cb -> [ [id: 'reference'], cb ] }
            .ifEmpty( [ [id: 'null'], [] ] )

        ARRIBA_VISUALISATION(
            ch_vis_input,
            ch_genes_gtf,
            ch_prot_vis,
            ch_cyto_vis
        )

    emit:
        fusions      = FUSION_ARRIBA.out.fusions        // [ meta, fusions.tsv ]
        fusions_fail = FUSION_ARRIBA.out.fusions_fail   // [ meta, fusions.discarded.tsv ]
        bam          = SAMTOOLS_SORT_ARRIBA.out.bam    // [ meta, sorted.bam, bai ]
        pdf          = ARRIBA_VISUALISATION.out.pdf    // [ meta, fusions.pdf ]
        star_log     = FUSION_ARRIBA.out.log_final      // [ meta, Log.final.out ]
}
