include { RENAME_FASTA_CONTIGS as RENAME_FASTA_CONTIGS_REF } from "../../../modules/local/rename_fasta_contigs/main.nf"
include { RENAME_DELIM_CONTIGS } from "../../../modules/local/rename_delim_contigs/main.nf"
include { WRITE_GENOME_CONFIG  } from "../../../modules/local/write_genome_config/main.nf"

include { STAR_GENOMEGENERATE } from "../../../modules/nf-core/star/genomegenerate/main.nf"

// Custom-genome-build helpers -- port of RENEE (classic) workflow/rules/build.smk's
// remaining rules, so a custom `--build` produces the same reference files a
// predefined conf/genomes/*.config entry ships with (see WRITE_GENOME_CONFIG).
include { BUILD_ANNOTATE       } from "../../../modules/local/build_annotate/main.nf"
include { BUILD_RRNA_LIST      } from "../../../modules/local/build_rrna_list/main.nf"
include { BUILD_TIN_REF        } from "../../../modules/local/build_tin_ref/main.nf"
include { BUILD_QUALIMAP_INFO  } from "../../../modules/local/build_qualimap_info/main.nf"
include { BUILD_RSEM_REF       } from "../../../modules/local/build_rsem_ref/main.nf"


workflow PREPARE_GENOME {

    take:
        // Directory to search for auto-detected Arriba fusion-calling reference
        // files during a custom genome build (may emit null); ignored for a
        // preset --genome, which sources its own fusion_* paths from conf/genomes.
        // Built in main.nf so a same-invocation --download_shared_resources run
        // can defer this until ARRIBA_DOWNLOAD actually finishes -- see the
        // .map{} usage below.
        ch_arriba_db_dir

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

            // Utils.optionalPathParam(g, key) will create a value channel with a value from the dict g using the key
            // or an empty channel when no key present
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

            ch_fasta_file = ch_fasta.map { meta, fa -> fa }
            ch_gtf_file   = ch_gtf.map   { meta, gtf -> gtf }

            def genome_name = params.genome ?: 'custom_genome'

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

            ch_organism = Channel.value( params.organism ?: 'custom' )

            // annotate.genes.txt / annotate.isoforms.txt / refFlat.txt / genes.ref.bed /
            // geneinfo.bed / karyobeds/ / karyoplot_gene_coordinates.txt all come
            // from one process (make_refFlat.py and make_geneinfo.py read each
            // other's output from the CWD with no arguments -- see
            // modules/local/build_annotate, which also folds in RENEE classic's
            // karyo_beds/karyo_coord rules since they too take only the GTF).
            // gene_info has no override param, so this always has to run; the
            // rest individually still defer to an explicit --annotate/--refflat/
            // etc. when one is given.
            BUILD_ANNOTATE(ch_gtf_file)
            ch_gene_info         = BUILD_ANNOTATE.out.gene_info
            ch_annotate          = params.annotate          ? Utils.optionalPathParam(params, 'annotate')          : BUILD_ANNOTATE.out.genes
            ch_annotate_isoforms = params.annotate_isoforms ? Utils.optionalPathParam(params, 'annotate_isoforms') : BUILD_ANNOTATE.out.isoforms
            ch_refflat           = params.refflat           ? Utils.optionalPathParam(params, 'refflat')           : BUILD_ANNOTATE.out.refflat
            ch_bed_ref           = params.bed_ref           ? Utils.optionalPathParam(params, 'bed_ref')           : BUILD_ANNOTATE.out.bed_ref
            ch_karyobeds         = params.karyobeds         ? Utils.optionalPathParam(params, 'karyobeds')         : BUILD_ANNOTATE.out.karyobeds
            ch_karyoploter       = params.karyoploter       ? Utils.optionalPathParam(params, 'karyoploter')       : BUILD_ANNOTATE.out.karyoploter

            if (params.qualimap_info) {
                ch_qualimap_info = Utils.optionalPathParam(params, 'qualimap_info')
            } else {
                ch_qualimap_info = BUILD_QUALIMAP_INFO(ch_fasta_file, ch_gtf_file).qualimap_info
            }

            if (params.rrna_list) {
                ch_rrna_list = Utils.optionalPathParam(params, 'rrna_list')
            } else {
                ch_rrna_list = BUILD_RRNA_LIST(ch_fasta_file, ch_gtf_file, genome_name).rrna_list
            }

            if (params.tin_ref) {
                ch_tin_ref = Utils.optionalPathParam(params, 'tin_ref')
            } else {
                ch_tin_ref = BUILD_TIN_REF(ch_gtf_file).tin_ref
            }

            // rsem_ref is a --reference *prefix* string, not a single file/dir (see
            // WRITE_GENOME_CONFIG). An explicit --rsem_ref points at an
            // already-built external reference, so BUILD_RSEM_REF only runs when
            // one isn't given.
            //
            // ch_rsem_ref (the subworkflow's own emitted value, used if this same
            // run goes on to RSEM quantification) always resolves to a real,
            // usable prefix. ch_rsem_ref_explicit is different on purpose: it's
            // what WRITE_GENOME_CONFIG uses to decide whether to stage the built
            // reference into the published genome bundle, so it must stay empty
            // in the built case -- an already-non-empty ch_rsem_ref there would
            // make WRITE_GENOME_CONFIG treat the just-built reference as an
            // external one and skip staging it, leaving the generated config
            // pointing at Nextflow's (ephemeral) work directory.
            ch_rsem_ref_explicit = Channel.value(params.rsem_ref ?: '')
            if (params.rsem_ref) {
                ch_rsem_ref     = Channel.value(params.rsem_ref)
                ch_rsem_ref_dir = Channel.value([])
            } else {
                ch_rsem_ref_dir = BUILD_RSEM_REF(ch_fasta_file, ch_gtf_file, genome_name).rsem_ref
                ch_rsem_ref     = ch_rsem_ref_dir.map { dir -> "${dir}/${genome_name}" }
            }
            // Arriba fusion-calling references. If not explicitly set, fall back to
            // auto-detecting them in ch_arriba_db_dir by matching params.genome
            // against known assembly names (hg19/hg38/mm10/mm39) -- this is how
            // `--build --shared_resources` output gets picked up for a custom build.
            // The glob-matching runs inside .map{} so it only fires once
            // ch_arriba_db_dir actually emits a value: for a same-invocation
            // `--download_shared_resources` run that's after ARRIBA_DOWNLOAD
            // finishes (see main.nf), not eagerly before it has run.
            ch_fusion_blacklist = ch_arriba_db_dir.map { dir ->
                def path = params.fusion_blacklist ?: Utils.arribaReferenceFile( dir, params.genome, 'blacklist_', '.tsv.gz' )
                path ? file( path, checkIfExists: true ) : []
            }
            ch_fusion_cytoband = ch_arriba_db_dir.map { dir ->
                def path = params.fusion_cytoband ?: Utils.arribaReferenceFile( dir, params.genome, 'cytobands_', '.tsv' )
                path ? file( path, checkIfExists: true ) : []
            }
            ch_fusion_protdomain = ch_arriba_db_dir.map { dir ->
                def path = params.fusion_protdomain ?: Utils.arribaReferenceFile( dir, params.genome, 'protein_domains_', '.gff3' )
                path ? file( path, checkIfExists: true ) : []
            }
            ch_fusion_known_fusions = ch_arriba_db_dir.map { dir ->
                def path = params.fusion_known_fusions ?: Utils.arribaReferenceFile( dir, params.genome, 'known_fusions_', '.tsv.gz' )
                path ? file( path, checkIfExists: true ) : []
            }

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
                ch_rsem_ref_explicit,
                ch_rsem_ref_dir,
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
