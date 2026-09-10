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
            ch_annotate          = Utils.optionalPathParam( g, 'annotate' )
            ch_annotate_isoforms = Utils.optionalPathParam( g, 'annotate_isoforms' )
            ch_refflat           = Utils.optionalPathParam( g, 'refflat' )
            ch_bed_ref           = Utils.optionalPathParam( g, 'bed_ref' )
            ch_qualimap_info     = Utils.optionalPathParam( g, 'qualimap_info' )
            ch_karyobeds         = Utils.optionalPathParam( g, 'karyobeds' )
            ch_karyoploter       = Utils.optionalPathParam( g, 'karyoploter' )
            ch_rsem_ref          = Channel.value( g.rsem_ref )
            ch_rrna_list         = Utils.optionalPathParam( g, 'rrna_list' )
            ch_tin_ref           = Utils.optionalPathParam( g, 'tin_ref' )

            // arriba vars
            ch_fusion_blacklist  = Utils.optionalPathParam( g, 'fusion_blacklist' )
            ch_fusion_cytoband   = Utils.optionalPathParam( g, 'fusion_cytoband' )
            ch_fusion_protdomain = Utils.optionalPathParam( g, 'fusion_protdomain' )

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
            ch_annotate          = Utils.optionalPathParam( params, 'annotate' )
            ch_annotate_isoforms = Utils.optionalPathParam( params, 'annotate_isoforms' )
            ch_refflat           = Utils.optionalPathParam( params, 'refflat' )
            ch_bed_ref           = Utils.optionalPathParam( params, 'bed_ref' )
            ch_qualimap_info     = Utils.optionalPathParam( params, 'qualimap_info' )
            ch_karyobeds         = Utils.optionalPathParam( params, 'karyobeds' )
            ch_karyoploter       = Utils.optionalPathParam( params, 'karyoploter' )
            ch_rsem_ref          = Channel.value( params.rsem_ref ?: '' )
            ch_rrna_list         = Utils.optionalPathParam( params, 'rrna_list' )
            ch_tin_ref           = Utils.optionalPathParam( params, 'tin_ref' )
            // Arriba fusion-calling references. If not explicitly set, fall back to
            // auto-detecting them in params.arriba_db_dir by matching params.genome
            // against known assembly names (hg19/hg38/mm10/mm39) -- this is how
            // `--build --shared_resources` output gets picked up for a custom build.
            ch_fusion_blacklist     = Utils.resolveOptionalPathParam( params, 'fusion_blacklist',
                Utils.arribaReferenceFile( params.arriba_db_dir, params.genome, 'blacklist_', '.tsv.gz' ) )
            ch_fusion_cytoband      = Utils.resolveOptionalPathParam( params, 'fusion_cytoband',
                Utils.arribaReferenceFile( params.arriba_db_dir, params.genome, 'cytobands_', '.tsv' ) )
            ch_fusion_protdomain    = Utils.resolveOptionalPathParam( params, 'fusion_protdomain',
                Utils.arribaReferenceFile( params.arriba_db_dir, params.genome, 'protein_domains_', '.gff3' ) )
            ch_fusion_known_fusions = Utils.resolveOptionalPathParam( params, 'fusion_known_fusions',
                Utils.arribaReferenceFile( params.arriba_db_dir, params.genome, 'known_fusions_', '.tsv.gz' ) )

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
