nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { PRESEQ_CCURVE }     from './modules/nf-core/preseq/ccurve/main'
include { HANDLE_PRESEQ_ERROR } from './modules/local/preseq/helperfunctions/main'
include { PARSE_PRESEQ_LOG } from './modules/local/preseq/helperfunctions/main'

// Modules
include { QUALIMAP_BAMQC } from './modules/nf-core/qualimap/bamqc/main'
include { SAMTOOLS_FLAGSTAT } from './modules/CCBR/samtools/flagstat/main.nf'
include { MULTIQC } from './modules/nf-core/multiqc/main'

// Subworkflows
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

        // build genome first, if set to build mode stop after this step
        PREPARE_GENOME()

        if (params.build_genome) {
            log.info "Build genome only mode enabled. Stopping workflow after genome preparation."
            prepare_genome_conf = PREPARE_GENOME.out.conf



        } else {
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

            // Estimate library complexity from mark-duplicated BAM (matches snakemake preseq rule)
// PICARD_INITIAL_QC.out.bam.view()
            PRESEQ_CCURVE(PICARD_INITIAL_QC.out.bam)

            // when preseq fails, write NAs for the stats that are calculated from its log
            PRESEQ_CCURVE.out.log
                .join(PICARD_INITIAL_QC.out.bam, remainder: true)
                .branch { meta, preseq_log, bam_tuple ->
                failed: preseq_log == null
                    return (tuple(meta, "nopresqlog"))
                succeeded: true
                    return (tuple(meta, preseq_log))
                }.set{ preseq_logs }
            preseq_logs.failed | HANDLE_PRESEQ_ERROR
            preseq_logs.succeeded | PARSE_PRESEQ_LOG
            PARSE_PRESEQ_LOG.out.nrf
                .concat(HANDLE_PRESEQ_ERROR.out.nrf)
                .set{ preseq_nrf }
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

            // MultiQC ------------------------------------------------------------------------
            ch_multiqc_files = Channel.empty()
                .mix(INITIAL_QC.out.fastqc_raw.map           { meta, files -> files }.flatten())
                .mix(INITIAL_QC.out.cutadapt_log.map         { meta, log   -> log   })
                .mix(INITIAL_QC.out.fastqc_trimmed.map       { meta, files -> files }.flatten())
                .mix(INITIAL_QC.out.fqscreen_1_txt.map       { meta, txt   -> txt   })
                .mix(INITIAL_QC.out.fqscreen_2_txt.map       { meta, txt   -> txt   })
                .mix(STAR_ALIGN.out.pass1_log.map            { meta, log   -> log   })
                .mix(STAR_ALIGN.out.pass2_log.map            { meta, log   -> log   })
                .mix(PICARD_COLLECTRNASEQMETRICS.out.metrics.map { meta, file -> file })
                .mix(QUALIMAP_BAMQC.out.results.map          { meta, dir   -> dir   })
                .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map      { meta, file  -> file  })
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
        prepare_genome_conf = params.build_genome ? prepare_genome_conf : Channel.empty()

        fastqc_raw     = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqc_raw
        fastqvalidator = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqvalidator
        cutadapt_reads = params.build_genome ? Channel.empty() : INITIAL_QC.out.cutadapt_reads
        cutadapt_log   = params.build_genome ? Channel.empty() : INITIAL_QC.out.cutadapt_log
        fastqc_trimmed = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqc_trimmed
        bbtools_ihist  = params.build_genome ? Channel.empty() : INITIAL_QC.out.bbtools_ihist
        fqscreen_1_txt = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_1_txt
        fqscreen_1_png = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_1_png
        fqscreen_2_txt = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_2_txt
        fqscreen_2_png = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_2_png

        star_pass1_sj             = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass1_sj
        star_pass1_log            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass1_log
        star_sjdb                 = params.build_genome ? Channel.empty() : STAR_ALIGN.out.sjdb
        star_pass2_log            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_log
        star_pass2_sj             = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_sj
        star_pass2_reads_per_gene = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_reads_per_gene
        star_pass2_bam            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_bam
        star_pass2_transcript_bam = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_transcript_bam

        picard_bam                = params.build_genome ? Channel.empty() : PICARD_INITIAL_QC.out.bam
        picard_bai                = params.build_genome ? Channel.empty() : PICARD_INITIAL_QC.out.bai
        picard_rna_metrics        = params.build_genome ? Channel.empty() : PICARD_COLLECTRNASEQMETRICS.out.metrics

        preseq_ccurve             = params.build_genome ? Channel.empty() : PRESEQ_CCURVE.out.c_curve
        preseq_log                = params.build_genome ? Channel.empty() : PRESEQ_CCURVE.out.log

        qualimap_results          = params.build_genome ? Channel.empty() : QUALIMAP_BAMQC.out.results

        flagstat                  = params.build_genome ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.flagstat
        flagstat_versions         = params.build_genome ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.versions

        arriba_fusions      = params.build_genome ? Channel.empty() : ARRIBA.out.fusions
        arriba_fusions_fail = params.build_genome ? Channel.empty() : ARRIBA.out.fusions_fail
        arriba_bam          = params.build_genome ? Channel.empty() : ARRIBA.out.bam
        arriba_pdf          = params.build_genome ? Channel.empty() : ARRIBA.out.pdf
        arriba_star_log     = params.build_genome ? Channel.empty() : ARRIBA.out.star_log

        multiqc_report      = params.build_genome ? Channel.empty() : MULTIQC.out.report
        multiqc_data        = params.build_genome ? Channel.empty() : MULTIQC.out.data
        preseq_nrf          = params.build_genome ? Channel.empty() : preseq_nrf
}

output {
    prepare_genome_conf { path { file -> "genome/" } }
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
    preseq_nrf          { path { meta, file -> 'preseq/' } }
}
