nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Modules
include { QUALIMAP_BAMQC } from './modules/nf-core/qualimap/bamqc/main'
include { SAMTOOLS_FLAGSTAT } from './modules/CCBR/samtools/flagstat/main.nf'

// Subworkflows
include { DOWNLOAD_DATABASES } from './subworkflows/local/download_databases/main.nf'
include { STAR_ALIGN }         from './subworkflows/local/star_align/main'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome/main.nf'
include { INITIAL_QC }         from './subworkflows/local/initial_qc/main'
include { CHECK_INPUT }     from './subworkflows/local/read_samples/main'
include { PICARD_INITIAL_QC }          from './subworkflows/local/picard_initial_qc/main'
include { PICARD_COLLECTRNASEQMETRICS } from './modules/nf-core/picard/collectrnaseqmetrics/main.nf'
include { arriba as ARRIBA }            from './subworkflows/local/arriba/main'





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

        // Download shared databases in either full build mode or resources-only mode.
        if (params.build_shared_resources_only || (params.build && params.shared_resources)) {
            log.info "Shared resources build enabled. Downloading FastQ Screen and Kraken2 databases."
            DOWNLOAD_DATABASES()
        }

        if (params.build_shared_resources_only) {
            log.info "Shared resources only mode enabled. Skipping genome preparation and sample analysis."
            prepare_genome_conf = Channel.empty()
        } else {
            // Prepare the genome for both genome-build and sample-analysis modes.
            PREPARE_GENOME()
        }

        if (params.build && !params.build_shared_resources_only) {
            log.info "Build genome only mode enabled. Stopping workflow after genome preparation."
            prepare_genome_conf = PREPARE_GENOME.out.conf
        } else if (!params.build_shared_resources_only) {
            log.info "Genome preparation complete. Continuing with workflow."
            prepare_genome_conf = Channel.empty() // dont save genome copy if not in build mode

            // Read samplesheet and emit channel of reads
            CHECK_INPUT(params.input)

            // Initial QC, trimming, and FastQ Screen steps
            INITIAL_QC(CHECK_INPUT.out.reads)


            // STAR alignment steps ----------------------------------------------------------

            STAR_ALIGN(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf
            )

            // post-alignment steps ----------------------------------------------------------
            PICARD_INITIAL_QC(STAR_ALIGN.out.pass2_bam)

            ch_fasta_path = PREPARE_GENOME.out.fasta.map { meta, fasta -> fasta }

            PICARD_COLLECTRNASEQMETRICS(
                PICARD_INITIAL_QC.out.bam,
                PREPARE_GENOME.out.refflat,
                ch_fasta_path,
                PREPARE_GENOME.out.rrna_list.ifEmpty([])
            )
            // QualiMap BAM QC ---------------------------------------------------------------
            ch_gtf_path = PREPARE_GENOME.out.genes_gtf.map { meta, gtf -> gtf }
            QUALIMAP_BAMQC(PICARD_INITIAL_QC.out.bam, ch_gtf_path)


            // SAMTOOLS_FLAGSTAT expects a tuple: [ meta, bam, bai ]
            picard_bam_bai_ch = PICARD_INITIAL_QC.out.bam.join(PICARD_INITIAL_QC.out.bai)

            SAMTOOLS_FLAGSTAT(picard_bam_bai_ch)


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

        fastq_screen_databases = (params.build_shared_resources_only || (params.shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.fastq_screen_databases : Channel.empty()
        kraken_databases       = (params.build_shared_resources_only || (params.shared_resources && params.build)) ? DOWNLOAD_DATABASES.out.kraken_databases       : Channel.empty()

        fastqc_raw     = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fastqc_raw
        fastqvalidator = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fastqvalidator
        cutadapt_reads = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.cutadapt_reads
        cutadapt_log   = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.cutadapt_log
        fastqc_trimmed = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fastqc_trimmed
        bbtools_ihist  = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.bbtools_ihist
        fqscreen_1_txt = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fqscreen_1_txt
        fqscreen_1_png = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fqscreen_1_png
        fqscreen_2_txt = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fqscreen_2_txt
        fqscreen_2_png = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.fqscreen_2_png
        kraken2_report                      = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.kraken2_report
        kraken2_classified_reads_assignment = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.kraken2_classified_reads_assignment
        kraken2_krona_html                  = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.kraken2_krona_html
        kraken2_db_dir                      = (params.build || params.build_shared_resources_only) ? Channel.empty() : INITIAL_QC.out.kraken2_db_dir

        star_pass1_sj             = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass1_sj
        star_pass1_log            = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass1_log
        star_sjdb                 = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.sjdb
        star_pass2_log            = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass2_log
        star_pass2_sj             = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass2_sj
        star_pass2_reads_per_gene = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass2_reads_per_gene
        star_pass2_bam            = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass2_bam
        star_pass2_transcript_bam = (params.build || params.build_shared_resources_only) ? Channel.empty() : STAR_ALIGN.out.pass2_transcript_bam

        picard_bam                = (params.build || params.build_shared_resources_only) ? Channel.empty() : PICARD_INITIAL_QC.out.bam
        picard_bai                = (params.build || params.build_shared_resources_only) ? Channel.empty() : PICARD_INITIAL_QC.out.bai

        qualimap_results          = (params.build || params.build_shared_resources_only) ? Channel.empty() : QUALIMAP_BAMQC.out.results

        flagstat                  = (params.build || params.build_shared_resources_only) ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.flagstat
        flagstat_versions         = (params.build || params.build_shared_resources_only) ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.versions

        arriba_fusions      = (params.build || params.build_shared_resources_only) ? Channel.empty() : ARRIBA.out.fusions
        arriba_fusions_fail = (params.build || params.build_shared_resources_only) ? Channel.empty() : ARRIBA.out.fusions_fail
        arriba_bam          = (params.build || params.build_shared_resources_only) ? Channel.empty() : ARRIBA.out.bam
        arriba_pdf          = (params.build || params.build_shared_resources_only) ? Channel.empty() : ARRIBA.out.pdf
        arriba_star_log     = (params.build || params.build_shared_resources_only) ? Channel.empty() : ARRIBA.out.star_log
}

output {
    // build outputs
    prepare_genome_conf {
        path { file -> "genome/" }
        mode 'copy'
        }
    fastq_screen_databases {
        path { meta, dir -> "${params.shared_resources}/fastq_screen_db/" }
        mode 'copy'
        }
    kraken_databases {
        path { meta, dir -> "${params.shared_resources}/" }
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

    picard_bam         { path { meta, file -> 'bams/' } }
    picard_bai         { path { meta, file -> 'bams/' } }
    picard_rna_metrics { path { meta, file -> 'picard/' } }

    qualimap_results { path { meta, dir -> "QualiMap/${meta.id}/" } }

    flagstat { path { meta, file -> 'log_files/' } }
    flagstat_versions { path { file -> "log_files/versions/${file.getParent().getFileName()}/" } }

    arriba_fusions      { path { meta, file -> 'fusions/' } }
    arriba_fusions_fail { path { meta, file -> 'fusions/' } }
    arriba_bam          { path { meta, bam, bai -> 'fusions/' } }
    arriba_pdf          { path { meta, file -> 'fusions/' } }
    arriba_star_log     { path { meta, file -> 'STAR_files/arriba/' } }
}
