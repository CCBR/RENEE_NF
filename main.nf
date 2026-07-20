nextflow.enable.dsl = 2

// Plugins
// Modules
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

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

            PICARD_COLLECTRNASEQMETRICS(
                PICARD_INITIAL_QC.out.bam,
                PREPARE_GENOME.out.refflat,
                [],
                PREPARE_GENOME.out.rrna_list
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

        arriba_fusions      = params.build_genome ? Channel.empty() : ARRIBA.out.fusions
        arriba_fusions_fail = params.build_genome ? Channel.empty() : ARRIBA.out.fusions_fail
        arriba_bam          = params.build_genome ? Channel.empty() : ARRIBA.out.bam
        arriba_pdf          = params.build_genome ? Channel.empty() : ARRIBA.out.pdf
        arriba_star_log     = params.build_genome ? Channel.empty() : ARRIBA.out.star_log
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

    arriba_fusions      { path { meta, file -> 'fusions/' } }
    arriba_fusions_fail { path { meta, file -> 'fusions/' } }
    arriba_bam          { path { meta, bam, bai -> 'fusions/' } }
    arriba_pdf          { path { meta, file -> 'fusions/' } }
    arriba_star_log     { path { meta, file -> 'STAR_files/arriba/' } }
}
