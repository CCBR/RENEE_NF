include { GTF2BED             } from "../../../modules/local/gtf2bed/main.nf"
include { RENAME_FASTA_CONTIGS as RENAME_FASTA_CONTIGS_REF } from "../../../modules/local/rename_fasta_contigs/main.nf"
include { RENAME_DELIM_CONTIGS } from "../../../modules/local/rename_delim_contigs/main.nf"
include { WRITE_GENOME_CONFIG  } from "../../../modules/local/write_genome_config/main.nf"

include { STAR_GENOMEGENERATE } from "../../../modules/nf-core/star/genomegenerate/main.nf"


workflow PREPARE_GENOME {

    main:

        ch_genome_conf = Channel.empty()

        if (params.genome && params.genomes && params.genomes[ params.genome ]) {

            def g = params.genomes[ params.genome ]

            ch_fasta     = Channel.value( [[id: params.genome], file( g.fasta, checkIfExists: true )] )
            ch_genes_gtf = Channel.value( [[id: params.genome], file( (params.star_gtf ?: g.genes_gtf), checkIfExists: true )] )
            ch_gene_info = Channel.value( file( g.gene_info, checkIfExists: true ) )

            def star_index_path = params.star_index ?: g.star_index
            if (star_index_path) {
                ch_star_index = Channel.value( [[id: params.genome],file( star_index_path, checkIfExists: true )] )
            } else {
                ch_star_index = STAR_GENOMEGENERATE(
                    ch_fasta,
                    ch_genes_gtf
                ).index
            }

            ch_organism          = Channel.value( g.organism )
            ch_annotate          = Channel.value( file( g.annotate,          checkIfExists: true ) )
            ch_annotate_isoforms = Channel.value( file( g.annotate_isoforms, checkIfExists: true ) )
            ch_refflat           = Channel.value( file( g.refflat,           checkIfExists: true ) )
            ch_bed_ref           = Channel.value( file( g.bed_ref,           checkIfExists: true ) )
            ch_qualimap_info     = Channel.value( file( g.qualimap_info,     checkIfExists: true ) )
            ch_karyobeds         = Channel.value( file( g.karyobeds,         checkIfExists: true ) )
            ch_karyoploter       = Channel.value( file( g.karyoploter,       checkIfExists: true ) )
            ch_rsem_ref          = Channel.value( g.rsem_ref )
            ch_rrna_list         = Channel.value( file( g.rrna_list,         checkIfExists: true ) )
            ch_tin_ref           = Channel.value( file( g.tin_ref,           checkIfExists: true ) )

            // arriba vars
            ch_fusion_blacklist     = g.fusion_blacklist     ? Channel.value( file( g.fusion_blacklist,     checkIfExists: true ) ) : Channel.empty()
            ch_fusion_cytoband      = g.fusion_cytoband      ? Channel.value( file( g.fusion_cytoband,      checkIfExists: true ) ) : Channel.empty()
            ch_fusion_protdomain    = g.fusion_protdomain    ? Channel.value( file( g.fusion_protdomain,    checkIfExists: true ) ) : Channel.empty()

            // fusion channel will often be empty, but will terminate the process if not a value channel
            ch_fusion_known_fusions = g.fusion_known_fusions ? Channel.value( file( g.fusion_known_fusions, checkIfExists: true ) ) : Channel.value([])

        } else if (params.genome_fasta && params.genes_gtf) {
            // If no genome config is provided, fall back to user-specified FASTA and GTF files.
            fasta_file = Channel.value( [[id: 'reference'], file( params.genome_fasta, checkIfExists: true )] )
            gtf_file   = Channel.value( [[id: 'reference'], file( params.genes_gtf, checkIfExists: true )] )

            if (params.rename_contigs) {

                contig_map = file( params.rename_contigs, checkIfExists: true )

                ch_fasta = RENAME_FASTA_CONTIGS_REF(
                    fasta_file,
                    contig_map
                ).fasta

                ch_gtf = RENAME_DELIM_CONTIGS(
                    gtf_file,
                    contig_map
                ).delim

            } else {

                ch_fasta = fasta_file
                ch_gtf = gtf_file
            }

            ch_genes_gtf = ch_gtf

            ch_gene_info = GTF2BED(ch_gtf).bed

            /*
             * Build STAR genome index for custom reference.
             *
             * You may need to adjust this call depending on the exact
             * input signature of your STAR_GENOMEGENERATE module.
             */
            ch_star_index = STAR_GENOMEGENERATE(
                ch_fasta,
                ch_genes_gtf
            ).index

            ch_organism          = Channel.value( params.organism ?: 'custom' )
            ch_annotate          = Channel.empty()
            ch_annotate_isoforms = Channel.empty()
            ch_refflat           = Channel.empty()
            ch_bed_ref           = Channel.empty()
            ch_qualimap_info     = Channel.empty()
            ch_karyobeds         = Channel.empty()
            ch_karyoploter       = Channel.empty()
            ch_rsem_ref          = Channel.value( params.rsem_ref ?: '' )
            ch_rrna_list         = Channel.empty()
            ch_tin_ref           = Channel.empty()
            ch_fusion_blacklist     = Channel.empty()
            ch_fusion_cytoband      = Channel.empty()
            ch_fusion_protdomain    = Channel.empty()
            ch_fusion_known_fusions = Channel.empty()

            WRITE_GENOME_CONFIG(
                ch_fasta,
                ch_genes_gtf,
                ch_gene_info,
                ch_star_index,
                ch_organism,
                ch_annotate.ifEmpty([]),
                ch_annotate_isoforms.ifEmpty([]),
                ch_refflat.ifEmpty([]),
                ch_bed_ref.ifEmpty([]),
                ch_qualimap_info.ifEmpty([]),
                ch_karyobeds.ifEmpty([]),
                ch_karyoploter.ifEmpty([]),
                ch_rsem_ref,
                ch_rrna_list.ifEmpty([]),
                ch_tin_ref.ifEmpty([]),
                ch_fusion_blacklist.ifEmpty([]),
                ch_fusion_cytoband.ifEmpty([]),
                ch_fusion_protdomain.ifEmpty([]),
                ch_fusion_known_fusions.ifEmpty([])
            )

            ch_genome_conf = WRITE_GENOME_CONFIG.out.conf.mix(
                WRITE_GENOME_CONFIG.out.files
            )

        } else {

            error """
            Either specify a genome in `conf/genomes.conf`,
            or specify both `--genome_fasta` and `--genes_gtf`
            to build a custom STAR reference for bulk RNA-seq.
            """
        }


    emit:
        // emits path unless otherwise specified
        fasta            = ch_fasta // tuple val(meta), path(fasta)
        genes_gtf        = ch_genes_gtf // tuple val(meta2), path(genes_gtf)
        gene_info        = ch_gene_info
        star_index       = ch_star_index // tuple val(meta4), path(star_index)
        organism         = ch_organism
        annotate         = ch_annotate
        annotate_isoforms = ch_annotate_isoforms
        refflat          = ch_refflat
        bed_ref          = ch_bed_ref
        qualimap_info    = ch_qualimap_info
        karyobeds        = ch_karyobeds
        karyoploter      = ch_karyoploter
        rsem_ref         = ch_rsem_ref
        rrna_list        = ch_rrna_list
        tin_ref          = ch_tin_ref
        fusion_blacklist     = ch_fusion_blacklist
        fusion_cytoband      = ch_fusion_cytoband
        fusion_protdomain    = ch_fusion_protdomain
        fusion_known_fusions = ch_fusion_known_fusions
        conf             = ch_genome_conf
}
