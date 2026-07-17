include { RSEM_CALCULATEEXPRESSION } from '../../../modules/nf-core/rsem/calculateexpression/main'
include { RSEM_MERGE               } from '../../../modules/local/rsem_merge/main'
include { RSEM_GENERATE_DATA_MATRIX } from '../../../modules/local/rsem_generate_data_matrix/main'

workflow RSEM {
    take:
        ch_transcript_bam  // channel: [ meta, bam ]  — toTranscriptome BAM from STAR pass2
        ch_strand_info     // channel: [ meta, path ] — strand.info from infer_experiment.py (or empty)
        ch_rsem_ref        // channel: val(rsem_ref_prefix string)
        ch_annotate        // channel: path(annotate.genes.txt)

    main:

        // Derive the RSEM index directory (parent of the ref prefix)
        // e.g. "/path/to/rsemref/hg38_36" → stage "/path/to/rsemref/"
        ch_rsem_index = ch_rsem_ref.map { rsem_ref ->
            file(rsem_ref).parent
        }

        // Pair transcript BAM with strand_info (by meta); fall back to null if empty
        ch_bam_with_strand = ch_transcript_bam
            .join(ch_strand_info, by: [0], remainder: true)

        // Compute meta.strandedness from strand.info content and build input channel
        ch_rsem_input = ch_bam_with_strand.map { meta, bam, strand_info ->
            def strandedness = 'none'
            if (strand_info) {
                def lines = strand_info.readLines().findAll { it.trim() }
                if (lines) {
                    def lastVal = lines.last().tokenize().last() as float
                    strandedness = lastVal > 0.75 ? 'reverse' : (lastVal < 0.25 ? 'forward' : 'none')
                }
            }
            [ meta + [strandedness: strandedness], bam ]
        }

        // Run RSEM calculate-expression (--alignments mode, BAM input)
        RSEM_CALCULATEEXPRESSION(
            ch_rsem_input,
            ch_rsem_index
        )

        // Collect all per-sample results for merging
        ch_genes_collected    = RSEM_CALCULATEEXPRESSION.out.counts_gene
            .map { meta, f -> f }
            .collect()

        ch_isoforms_collected = RSEM_CALCULATEEXPRESSION.out.counts_transcript
            .map { meta, f -> f }
            .collect()

        // Merge per-sample results into count matrices using RENEE's Python script
        RSEM_MERGE(
            ch_genes_collected,
            ch_isoforms_collected,
            ch_annotate,
            '.'
        )

        // Generate raw data matrices using rsem-generate-data-matrix
        RSEM_GENERATE_DATA_MATRIX(
            ch_genes_collected,
            ch_isoforms_collected
        )

    emit:
        genes_results    = RSEM_CALCULATEEXPRESSION.out.counts_gene
        isoforms_results = RSEM_CALCULATEEXPRESSION.out.counts_transcript
        stat             = RSEM_CALCULATEEXPRESSION.out.stat
        gene_counts      = RSEM_MERGE.out.gene_counts
        gene_fpkm        = RSEM_MERGE.out.gene_fpkm
        gene_tpm         = RSEM_MERGE.out.gene_tpm
        isoform_counts   = RSEM_MERGE.out.isoform_counts
        isoform_fpkm     = RSEM_MERGE.out.isoform_fpkm
        isoform_tpm      = RSEM_MERGE.out.isoform_tpm
        reformatted      = RSEM_MERGE.out.reformatted
        gene_matrix      = RSEM_GENERATE_DATA_MATRIX.out.gene_matrix
        isoform_matrix   = RSEM_GENERATE_DATA_MATRIX.out.isoform_matrix
}
