nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


// Modules
include { MULTIQC }                     from './modules/nf-core/multiqc/main'
include { BAM2STRANDEDBW }              from './modules/local/bam2strandedbw/main'
include { FC_LANE }                     from './modules/local/fc_lane/main'
include { MULTIQCPARSER }               from './modules/local/multiqcparser/main'
include { RNA_REPORT }                  from './modules/local/rna_report/main'

// Subworkflows
include { DOWNLOAD_DATABASES }  from './subworkflows/local/download_databases/main.nf'
include { STAR_ALIGN_WORKFLOW } from './subworkflows/local/star_align/main'
include { PREPARE_GENOME }      from './subworkflows/local/prepare_genome/main.nf'
include { INITIAL_QC }          from './subworkflows/local/initial_qc/main'
include { CHECK_INPUT }         from './subworkflows/local/read_samples/main'
include { PICARD_INITIAL_QC }   from './subworkflows/local/picard_initial_qc/main'
include { POST_ALIGNMENT_QC }   from './subworkflows/local/post_alignment_qc/main'
include { RSEM }                from './subworkflows/local/rsem/main'
include { arriba as ARRIBA }   from './subworkflows/local/arriba/main'


workflow version {
    println "RENEE_NF ${workflow.manifest.version}"
}

workflow LOG {
    log.info """\
            RENEE_NF $workflow.manifest.version
            =============
            cmd line     : $workflow.commandLine
            start time   : $workflow.start
            launchDir    : $workflow.launchDir
            input        : ${params.input}
            genome       : ${params.genome}
            """
            .stripIndent()
    log.info paramsSummaryLog(workflow)
}

workflow {
    main:
        LOG()
        validateParameters()

        if (params.build && params.build_shared_resources_only) {
            error "Parameters --build and --build_shared_resources_only are mutually exclusive."
        }

        analysis_mode = !params.build && !params.build_shared_resources_only

        // Download shared databases in either full build mode or resources-only mode.
        // download_shared_resources is the on/off switch; the destination is always
        // <outputDir>/shared_resources (Utils.sharedResourcesDir), which INITIAL_QC
        // also auto-detects from for an analysis run reusing the same --outputDir.
        if (params.build_shared_resources_only || (params.build && params.download_shared_resources)) {
            log.info "Shared resources build enabled. Downloading FastQ Screen and Kraken2 databases."
            DOWNLOAD_DATABASES()
        }

        // Directory to search for auto-detected Arriba fusion-calling reference
        // files during a custom genome build. An explicit --arriba_db_dir always
        // wins (e.g. a directory populated by a prior separate run); otherwise,
        // when this same invocation also downloads shared resources, defer until
        // ARRIBA_DOWNLOAD actually finishes -- .collect() only emits once its
        // source channel closes -- rather than letting PREPARE_GENOME's file-glob
        // lookup run before the download has happened.
        if (params.arriba_db_dir) {
            ch_arriba_db_dir = Channel.value(params.arriba_db_dir)
        } else if (params.build_shared_resources_only || (params.build && params.download_shared_resources)) {
            ch_arriba_db_dir = DOWNLOAD_DATABASES.out.arriba_database
                .collect()
                .map { files ->
                    def first_file = files.flatten().find()
                    first_file ? first_file.parent.toString() : null
                }
        } else {
            ch_arriba_db_dir = Channel.value([]) // Use ([]) instead of (null)
        }

        if (params.build_shared_resources_only) {
            log.info "Shared resources only mode enabled. Skipping genome preparation and sample analysis."
            prepare_genome_conf = Channel.empty()
        } else {
            // Prepare the genome for both genome-build and sample-analysis modes.
            PREPARE_GENOME(ch_arriba_db_dir)
        }

        if (params.build && !params.build_shared_resources_only) {
            log.info "Build genome only mode enabled. Stopping workflow after genome preparation."
            prepare_genome_conf = PREPARE_GENOME.out.conf
        } else if (analysis_mode) {
            log.info "Genome preparation complete. Continuing with workflow."
            prepare_genome_conf = Channel.empty() // dont save genome copy if not in build mode

            // Read samplesheet and emit channel of reads
            CHECK_INPUT(params.input)

            // Initial QC, trimming, and FastQ Screen steps
            INITIAL_QC(CHECK_INPUT.out.reads)


            // STAR alignment steps ----------------------------------------------------------

            STAR_ALIGN_WORKFLOW(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf
            )

            // post-alignment steps ----------------------------------------------------------
            PICARD_INITIAL_QC(STAR_ALIGN_WORKFLOW.out.pass2_bam)

            POST_ALIGNMENT_QC(
                PICARD_INITIAL_QC.out.bam,
                PICARD_INITIAL_QC.out.bai,
                PREPARE_GENOME.out.bed_ref,
                PREPARE_GENOME.out.tin_ref,
                PREPARE_GENOME.out.refflat,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.rrna_list,
                PREPARE_GENOME.out.genes_gtf
            )

            // RSEM quantification -----------------------------------------------
            RSEM(
                STAR_ALIGN_WORKFLOW.out.pass2_transcript_bam,
                POST_ALIGNMENT_QC.out.infer_experiment,
                PREPARE_GENOME.out.rsem_ref,
                PREPARE_GENOME.out.annotate
            )
            // BAM to stranded BigWig files
            BAM2STRANDEDBW(
                PICARD_INITIAL_QC.out.bam
                    .join(PICARD_INITIAL_QC.out.bai)
                    .join(POST_ALIGNMENT_QC.out.infer_experiment)
            )


            // Arriba gene-fusion calling (only when genome supplies a blacklist) ----------

            ARRIBA(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fusion_blacklist,
                PREPARE_GENOME.out.fusion_known_fusions,
                PREPARE_GENOME.out.fusion_cytoband,
                PREPARE_GENOME.out.fusion_protdomain
            )

            // MultiQC ------------------------------------------------------------------------
            ch_multiqc_files = Channel.empty()
                .mix(INITIAL_QC.out.fastqc_raw.map           { meta, files -> files }.flatten())
                .mix(INITIAL_QC.out.cutadapt_log.map         { meta, log   -> log   })
                .mix(INITIAL_QC.out.fastqc_trimmed.map       { meta, files -> files }.flatten())
                .mix(INITIAL_QC.out.fqscreen_1_txt.map       { meta, txt   -> txt   })
                .mix(INITIAL_QC.out.fqscreen_2_txt.map       { meta, txt   -> txt   })
                .mix(INITIAL_QC.out.kraken2_report.map       { meta, file  -> file  })
                .mix(STAR_ALIGN_WORKFLOW.out.pass1_log.map   { meta, log   -> log   })
                .mix(STAR_ALIGN_WORKFLOW.out.pass2_log.map   { meta, log   -> log   })
                .mix(PICARD_INITIAL_QC.out.metrics.map              { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.picard_rna_metrics.map  { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.qualimap_results.map    { meta, dir  -> dir  })
                .mix(POST_ALIGNMENT_QC.out.flagstat.map            { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.read_distribution.map   { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.infer_experiment.map    { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.inner_distance_freq.map { meta, file -> file })
                .mix(POST_ALIGNMENT_QC.out.tin_txt.map             { meta, file -> file })
                .mix(RSEM.out.genes_results.map              { meta, file  -> file  })
                .mix(POST_ALIGNMENT_QC.out.preseq_ccurve.map { meta, file  -> file  })
                .collect()

            // multiqc_config = channel.value(file('conf/multiqc_config.yaml'))
            MULTIQC(
                ch_multiqc_files.map { files -> [
                    [id: 'multiqc'], // meta
                    files, // files
                    file(params.multiQC_config), // config
                    file(params.multiQC_logo), //logo
                    [], // replace_names
                    []] //sample names TSV
                }
            )

            // Parse MultiQC outputs and render the RNA QC report ----------------------------
            FC_LANE(
                CHECK_INPUT.out.reads.map { meta, reads -> tuple(meta, reads[0]) }
            )

            ch_inner_distance_files = POST_ALIGNMENT_QC.out.inner_distance_freq
                .map { meta, file -> file }
                .collect()
                .map { files -> [files] }

            ch_tin_summary_files = POST_ALIGNMENT_QC.out.tin_txt
                .map { meta, file -> file }
                .collect()
                .map { files -> [files] }

            ch_fastq_info_files = FC_LANE.out.fqinfo
                .map { meta, file -> file }
                .collect()
                .map { files -> [files] }

            ch_multiqcparser_input = MULTIQC.out.data
                .combine(ch_inner_distance_files)
                .combine(ch_tin_summary_files)
                .combine(ch_fastq_info_files)
                .map { meta, data, inner_distance, tin_summary, fastq_info ->
                    tuple(meta, data, inner_distance, tin_summary, fastq_info)
                }

            MULTIQCPARSER(ch_multiqcparser_input)

            ch_tin_matrix_files = POST_ALIGNMENT_QC.out.tin_xls
                .map { meta, file -> file }
                .collect()
                .map { files -> [files] }

            ch_rna_report_input = MULTIQCPARSER.out.matrix
                .combine(RSEM.out.reformatted.map { counts -> [counts] })
                .combine(ch_tin_matrix_files)
                .map { meta, qc, counts, tins ->
                    tuple(
                        [id: 'rna_report'],
                        counts,
                        tins,
                        qc,
                        file("${projectDir}/assets/rNA_flowcells.Rmd")
                    )
                }

            RNA_REPORT(ch_rna_report_input)
        }
        workflow.onComplete = {
            if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
                def message = Utils.spooker(workflow)
                if (message) {
                    println message
                }
            }
        }

    publish:
        // In build genome mode, only publish the genome conf file, otherwise publish all outputs
        prepare_genome_conf = params.build ? prepare_genome_conf : Channel.empty()

        fastq_screen_databases = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.fastq_screen_databases : Channel.empty()
        kraken_databases       = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.kraken_databases       : Channel.empty()
        fastq_screen_conf1     = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.fastq_screen_conf1     : Channel.empty()
        fastq_screen_conf2     = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.fastq_screen_conf2     : Channel.empty()
        arriba_database        = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.arriba_database        : Channel.empty()
        shared_resources_conf  = (params.build_shared_resources_only || (params.download_shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.conf                  : Channel.empty()

        fastqc_raw     = analysis_mode ? INITIAL_QC.out.fastqc_raw : Channel.empty()
        fastqvalidator = analysis_mode ? INITIAL_QC.out.fastqvalidator : Channel.empty()
        cutadapt_reads = analysis_mode ? INITIAL_QC.out.cutadapt_reads : Channel.empty()
        cutadapt_log   = analysis_mode ? INITIAL_QC.out.cutadapt_log : Channel.empty()
        fastqc_trimmed = analysis_mode ? INITIAL_QC.out.fastqc_trimmed : Channel.empty()
        bbtools_ihist  = analysis_mode ? INITIAL_QC.out.bbtools_ihist : Channel.empty()
        fqscreen_1_txt = analysis_mode ? INITIAL_QC.out.fqscreen_1_txt : Channel.empty()
        fqscreen_1_png = analysis_mode ? INITIAL_QC.out.fqscreen_1_png : Channel.empty()
        fqscreen_2_txt = analysis_mode ? INITIAL_QC.out.fqscreen_2_txt : Channel.empty()
        fqscreen_2_png = analysis_mode ? INITIAL_QC.out.fqscreen_2_png : Channel.empty()
        kraken2_report                      = analysis_mode ? INITIAL_QC.out.kraken2_report : Channel.empty()
        kraken2_classified_reads_assignment = analysis_mode ? INITIAL_QC.out.kraken2_classified_reads_assignment : Channel.empty()
        kraken2_krona_html                  = analysis_mode ? INITIAL_QC.out.kraken2_krona_html : Channel.empty()
        kraken2_db_dir                      = analysis_mode ? INITIAL_QC.out.kraken2_db_dir : Channel.empty()

        star_pass1_sj             = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass1_sj : Channel.empty()
        star_pass1_log            = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass1_log : Channel.empty()
        star_sjdb                 = analysis_mode ? STAR_ALIGN_WORKFLOW.out.sjdb : Channel.empty()
        star_pass2_log            = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass2_log : Channel.empty()
        star_pass2_sj             = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass2_sj : Channel.empty()
        star_pass2_reads_per_gene = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass2_reads_per_gene : Channel.empty()
        star_pass2_bam            = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass2_bam : Channel.empty()
        star_pass2_transcript_bam = analysis_mode ? STAR_ALIGN_WORKFLOW.out.pass2_transcript_bam : Channel.empty()

        picard_bam             = analysis_mode ? PICARD_INITIAL_QC.out.bam : Channel.empty()
        picard_bai             = analysis_mode ? PICARD_INITIAL_QC.out.bai : Channel.empty()
        picard_markdup_metrics = analysis_mode ? PICARD_INITIAL_QC.out.metrics : Channel.empty()
        picard_rna_metrics     = analysis_mode ? POST_ALIGNMENT_QC.out.picard_rna_metrics : Channel.empty()

        qualimap_results = analysis_mode ? POST_ALIGNMENT_QC.out.qualimap_results : Channel.empty()

        flagstat          = analysis_mode ? POST_ALIGNMENT_QC.out.flagstat          : Channel.empty()
        flagstat_versions = analysis_mode ? POST_ALIGNMENT_QC.out.flagstat_versions : Channel.empty()

        arriba_fusions      = analysis_mode ? ARRIBA.out.fusions : Channel.empty()
        arriba_fusions_fail = analysis_mode ? ARRIBA.out.fusions_fail : Channel.empty()
        arriba_bam          = analysis_mode ? ARRIBA.out.bam : Channel.empty()
        arriba_pdf          = analysis_mode ? ARRIBA.out.pdf : Channel.empty()
        arriba_star_log     = analysis_mode ? ARRIBA.out.star_log : Channel.empty()

        rsem_genes_results    = analysis_mode ? RSEM.out.genes_results    : Channel.empty()
        rsem_isoforms_results = analysis_mode ? RSEM.out.isoforms_results : Channel.empty()
        rsem_gene_counts      = analysis_mode ? RSEM.out.gene_counts      : Channel.empty()
        rsem_gene_fpkm        = analysis_mode ? RSEM.out.gene_fpkm        : Channel.empty()
        rsem_gene_tpm         = analysis_mode ? RSEM.out.gene_tpm         : Channel.empty()
        rsem_isoform_counts   = analysis_mode ? RSEM.out.isoform_counts   : Channel.empty()
        rsem_isoform_fpkm     = analysis_mode ? RSEM.out.isoform_fpkm     : Channel.empty()
        rsem_isoform_tpm      = analysis_mode ? RSEM.out.isoform_tpm      : Channel.empty()
        rsem_reformatted      = analysis_mode ? RSEM.out.reformatted      : Channel.empty()
        rsem_gene_matrix      = analysis_mode ? RSEM.out.gene_matrix      : Channel.empty()
        rsem_isoform_matrix   = analysis_mode ? RSEM.out.isoform_matrix   : Channel.empty()

        rseqc_infer_experiment       = analysis_mode ? POST_ALIGNMENT_QC.out.infer_experiment       : Channel.empty()
        rseqc_read_distribution      = analysis_mode ? POST_ALIGNMENT_QC.out.read_distribution      : Channel.empty()
        rseqc_inner_distance_freq    = analysis_mode ? POST_ALIGNMENT_QC.out.inner_distance_freq    : Channel.empty()
        rseqc_inner_distance_dist    = analysis_mode ? POST_ALIGNMENT_QC.out.inner_distance_dist    : Channel.empty()
        rseqc_inner_distance_rscript = analysis_mode ? POST_ALIGNMENT_QC.out.inner_distance_rscript : Channel.empty()
        rseqc_tin_txt                = analysis_mode ? POST_ALIGNMENT_QC.out.tin_txt                : Channel.empty()
        rseqc_tin_xls                = analysis_mode ? POST_ALIGNMENT_QC.out.tin_xls                : Channel.empty()

        bam2bw_fwd = analysis_mode ? BAM2STRANDEDBW.out.fwd_bw : Channel.empty()
        bam2bw_rev = analysis_mode ? BAM2STRANDEDBW.out.rev_bw : Channel.empty()

        preseq_ccurve = analysis_mode ? POST_ALIGNMENT_QC.out.preseq_ccurve : Channel.empty()
        preseq_log    = analysis_mode ? POST_ALIGNMENT_QC.out.preseq_log    : Channel.empty()
        preseq_nrf    = analysis_mode ? POST_ALIGNMENT_QC.out.preseq_nrf    : Channel.empty()

        multiqc_report = analysis_mode ? MULTIQC.out.report : Channel.empty()
        multiqc_data   = analysis_mode ? MULTIQC.out.data   : Channel.empty()
        fastq_info     = analysis_mode ? FC_LANE.out.fqinfo  : Channel.empty()

        multiqc_matrix       = analysis_mode ? MULTIQCPARSER.out.matrix          : Channel.empty()
        rseqc_inner_distances = analysis_mode ? MULTIQCPARSER.out.inner_distances : Channel.empty()
        rseqc_median_tin     = analysis_mode ? MULTIQCPARSER.out.median_tin       : Channel.empty()
        fastq_flowcell_lanes = analysis_mode ? MULTIQCPARSER.out.flowcell_lanes   : Channel.empty()
        rna_report           = analysis_mode ? RNA_REPORT.out.html                : Channel.empty()

}

output {
    // build outputs
    prepare_genome_conf {
        path { file -> "genome/" }
        mode 'copy'
        }
    fastq_screen_databases {
        path { meta, dir -> "shared_resources/fastq_screen_db/" }
        mode 'copy'
        }
    kraken_databases {
        path { meta, dir -> "shared_resources/" }
        mode 'copy'
        }
    fastq_screen_conf1 {
        path { file -> "shared_resources/fastq_screen_db/" }
        mode 'copy'
        }
    fastq_screen_conf2 {
        path { file -> "shared_resources/fastq_screen_db/" }
        mode 'copy'
        }
    arriba_database {
        // Version pin lives in modules/nf-core/arriba/download/main.nf (not
        // exposed as a param there); keep this path in sync with it.
        // Each of the 4 mixed emits (blacklist/cytobands/protein_domains/
        // known_fusions) is a glob match across all genome builds, so this
        // closure receives a list of files, not a single one.
        path { files -> "shared_resources/arriba_v2.5.0/database/" }
        mode 'copy'
        }
    shared_resources_conf {
        // Sibling of shared_resources/, mirroring prepare_genome_conf's
        // <genome name>.config sitting alongside genome/<genome name>/.
        path { file -> "./" }
        mode 'copy'
        }

    // analysis outputs
    fastqc_raw { path { meta, file -> "fastqc/raw/" } }
    fastqvalidator { path { meta, file -> "fastqvalidator/${meta.id}/" } }
    cutadapt_reads { path { meta, reads -> "cutadapt/${meta.id}/" } }
    cutadapt_log { path { meta, log -> "cutadapt/${meta.id}/" } }
    fastqc_trimmed { path { meta, file -> "fastqc/trimmed/" } }
    bbtools_ihist { path { meta, ihist -> "bbtools/${meta.id}/" }}
    fqscreen_1_txt  { path { meta, file  -> "FQscreen/" } }
    fqscreen_1_png  { path { meta, file  -> "FQscreen/" } }
    fqscreen_2_txt  { path { meta, file  -> "FQscreen2/" } }
    fqscreen_2_png  { path { meta, file  -> "FQscreen2/" } }
    kraken2_report                      { path { meta, file -> "kraken2/" } }
    kraken2_classified_reads_assignment { path { meta, file -> "kraken2/" } }
    kraken2_krona_html                  { path { meta, file -> "kraken2/" } }
    kraken2_db_dir                      { path { dir -> "./" } }

    star_pass1_sj { path { meta, file -> 'STAR_files/pass1/' } }
    star_pass1_log { path { meta, file -> 'STAR_files/pass1/' } }
    star_sjdb { path { file -> 'STAR_files/pass1/' } }
    star_pass2_log { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_sj { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_reads_per_gene { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_bam { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_transcript_bam { path { meta, file -> 'bams/' } }

    rsem_genes_results    { path { meta, file -> 'DEG_ALL/' } }
    rsem_isoforms_results { path { meta, file -> 'DEG_ALL/' } }
    rsem_gene_counts      { path { file -> 'DEG_ALL/' } }
    rsem_gene_fpkm        { path { file -> 'DEG_ALL/' } }
    rsem_gene_tpm         { path { file -> 'DEG_ALL/' } }
    rsem_isoform_counts   { path { file -> 'DEG_ALL/' } }
    rsem_isoform_fpkm     { path { file -> 'DEG_ALL/' } }
    rsem_isoform_tpm      { path { file -> 'DEG_ALL/' } }
    rsem_reformatted      { path { file -> 'DEG_ALL/' } }
    rsem_gene_matrix      { path { file -> 'DEG_ALL/' } }
    rsem_isoform_matrix   { path { file -> 'DEG_ALL/' } }

    picard_bam             { path { meta, file -> 'bams/' } }
    picard_bai             { path { meta, file -> 'bams/' } }
    picard_markdup_metrics { path { meta, file -> 'picard/' } }
    picard_rna_metrics     { path { meta, file -> 'picard/' } }

    rseqc_infer_experiment       { path { meta, file -> 'RSeQC/' } }
    rseqc_read_distribution      { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_freq    { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_dist    { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_rscript { path { meta, file -> 'RSeQC/' } }
    rseqc_tin_txt                { path { meta, file -> 'RSeQC/' } }
    rseqc_tin_xls                { path { meta, file -> 'RSeQC/' } }

    bam2bw_fwd { path { meta, file -> 'bigwigs/' } }
    bam2bw_rev { path { meta, file -> 'bigwigs/' } }

    preseq_ccurve { path { meta, file -> 'preseq/' } }
    preseq_log    { path { meta, file -> 'preseq/' } }

    qualimap_results { path { meta, dir -> "QualiMap/${meta.id}/" } }

    flagstat { path { meta, file -> 'log_files/' } }
    flagstat_versions { path { file -> "log_files/versions/${file.getParent().getFileName()}/" } }

    arriba_fusions      { path { meta, file -> 'fusions/' } }
    arriba_fusions_fail { path { meta, file -> 'fusions/' } }
    arriba_bam          { path { meta, bam, bai -> 'fusions/' } }
    arriba_pdf          { path { meta, file -> 'fusions/' } }
    arriba_star_log     { path { meta, file -> 'STAR_files/arriba/' } }

    multiqc_report { path { meta, file -> 'Reports/' } }
    multiqc_data   { path { meta, dir  -> 'Reports/' } }
    fastq_info { path { meta, file -> 'rawQC/' } }
    multiqc_matrix { path { meta, file -> 'Reports/' } }
    rseqc_inner_distances { path { meta, file -> 'Reports/' } }
    rseqc_median_tin { path { meta, file -> 'Reports/' } }
    fastq_flowcell_lanes { path { meta, file -> 'Reports/' } }
    rna_report { path { meta, file -> 'Reports/' } }
    preseq_nrf          { path { meta, file -> 'preseq/' } }
}
