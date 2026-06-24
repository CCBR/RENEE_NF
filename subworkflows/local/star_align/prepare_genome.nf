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

            ch_fasta = Channel.fromPath(
                params.genomes[ params.genome ].fasta,
                checkIfExists: true
            )

            ch_genes_gtf = Channel.fromPath(
                params.genomes[ params.genome ].genes_gtf,
                checkIfExists: true
            )

            ch_gene_info = Channel.fromPath(
                params.genomes[ params.genome ].gene_info,
                checkIfExists: true
            )

            if (params.genomes[ params.genome ].star_index) {
                ch_star_index = Channel.fromPath(
                    params.genomes[ params.genome ].star_index,
                    type: 'dir',
                    checkIfExists: true
                )
            } else {
                ch_star_index = STAR_GENOMEGENERATE(
                    ch_fasta.map { fa -> [ [:], fa ] },
                    ch_genes_gtf.map { gtf -> [ [:], gtf ] }
                ).index.map { meta, idx -> idx }
            }

            ch_bioc_txdb = Channel.value(
                params.genomes[ params.genome ].bioc_txdb
            )

            ch_bioc_annot = Channel.value(
                params.genomes[ params.genome ].bioc_annot
            )

        } else if (params.genome_fasta && params.genes_gtf) {

            fasta_file = Channel.fromPath(
                params.genome_fasta,
                checkIfExists: true
            )

            gtf_file = file(
                params.genes_gtf,
                checkIfExists: true
            )

            if (params.rename_contigs) {

                contig_map = file(
                    params.rename_contigs,
                    checkIfExists: true
                )

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
                ch_gtf = Channel.fromPath(
                    params.genes_gtf,
                    checkIfExists: true
                )
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

            ch_bioc_txdb = Channel.value(params.bioc_txdb)
            ch_bioc_annot = Channel.value(params.bioc_annot)

            /*
             * Optional.
             * This call may need to be changed if WRITE_GENOME_CONFIG
             * still expects the old full set of genome-prep inputs.
             */
            WRITE_GENOME_CONFIG(
                ch_fasta,
                ch_genes_gtf,
                ch_gene_info,
                ch_star_index,
                ch_bioc_txdb,
                ch_bioc_annot
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

        fasta = ch_fasta
        genes_gtf  = ch_genes_gtf.map  { gtf -> [ [id: params.genome], gtf ] }
        gene_info  = ch_gene_info
        star_index = ch_star_index.map { idx -> [ [id: params.genome], idx ] }
        bioc_txdb  = ch_bioc_txdb
        bioc_annot = ch_bioc_annot
        conf       = ch_genome_conf
}
