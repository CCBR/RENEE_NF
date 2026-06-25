include {
    RENAME_FASTA_CONTIGS as RENAME_FASTA_CONTIGS_REF
    RENAME_DELIM_CONTIGS
    GTF2BED
    WRITE_GENOME_CONFIG
} from "../../../modules/local/prepare_genome.nf"

include { STAR_GENOMEGENERATE } from "../../../modules/nf-core/star/genomegenerate/main.nf"


workflow prepare_genome {

    main:

        ch_genome_conf = Channel.empty()

        if (params.genomes[ params.genome ]) {

            def g = params.genomes[ params.genome ]

            ch_fasta     = Channel.fromPath( g.fasta,     checkIfExists: true )
            ch_genes_gtf = Channel.fromPath( g.genes_gtf, checkIfExists: true )
            ch_gene_info = Channel.fromPath( g.gene_info, checkIfExists: true )

            if (g.star_index) {
                ch_star_index = Channel.fromPath( g.star_index, type: 'dir', checkIfExists: true )
            } else {
                ch_star_index = STAR_GENOMEGENERATE(
                    ch_fasta.map { fa -> [ [:], fa ] },
                    ch_genes_gtf.map { gtf -> [ [:], gtf ] }
                ).index.map { meta, idx -> idx }
            }

            ch_organism          = Channel.value( g.organism )
            ch_annotate          = Channel.fromPath( g.annotate,          checkIfExists: true )
            ch_annotate_isoforms = Channel.fromPath( g.annotate_isoforms, checkIfExists: true )
            ch_refflat           = Channel.fromPath( g.refflat,           checkIfExists: true )
            ch_bed_ref           = Channel.fromPath( g.bed_ref,           checkIfExists: true )
            ch_qualimap_info     = Channel.fromPath( g.qualimap_info,     checkIfExists: true )
            ch_karyobeds         = Channel.fromPath( g.karyobeds,         checkIfExists: true )
            ch_karyoploter       = Channel.fromPath( g.karyoploter,       checkIfExists: true )
            ch_rsem_ref          = Channel.value( g.rsem_ref )
            ch_rrna_list         = Channel.fromPath( g.rrna_list,         checkIfExists: true )
            ch_tin_ref           = Channel.fromPath( g.tin_ref,           checkIfExists: true )

            ch_fusion_blacklist  = g.fusion_blacklist  ? Channel.fromPath( g.fusion_blacklist,  checkIfExists: true ) : Channel.empty()
            ch_fusion_cytoband   = g.fusion_cytoband   ? Channel.fromPath( g.fusion_cytoband,   checkIfExists: true ) : Channel.empty()
            ch_fusion_protdomain = g.fusion_protdomain ? Channel.fromPath( g.fusion_protdomain, checkIfExists: true ) : Channel.empty()

        } else if (params.genome_fasta && params.genes_gtf) {

            fasta_file = Channel.fromPath( params.genome_fasta, checkIfExists: true )
            gtf_file   = file( params.genes_gtf, checkIfExists: true )

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
                ch_gtf = Channel.fromPath( params.genes_gtf, checkIfExists: true )
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
                ch_fasta.map { fa -> [ [:], fa ] },
                ch_genes_gtf.map { gtf -> [ [:], gtf ] }
            ).index.map { meta, idx -> idx }

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
            ch_fusion_blacklist  = Channel.empty()
            ch_fusion_cytoband   = Channel.empty()
            ch_fusion_protdomain = Channel.empty()

            WRITE_GENOME_CONFIG(
                ch_fasta,
                ch_genes_gtf,
                ch_gene_info,
                ch_star_index
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

        fasta            = ch_fasta
        genes_gtf        = ch_genes_gtf.map  { gtf -> [ [id: params.genome], gtf ] }
        gene_info        = ch_gene_info
        star_index       = ch_star_index.map { idx -> [ [id: params.genome], idx ] }
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
        fusion_blacklist  = ch_fusion_blacklist
        fusion_cytoband   = ch_fusion_cytoband
        fusion_protdomain = ch_fusion_protdomain
        conf             = ch_genome_conf
}
