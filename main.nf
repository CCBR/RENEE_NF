nextflow.enable.dsl = 2

// Plugins
// Modules
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Subworkflows
include { star_align_workflow as STAR_ALIGN }         from './subworkflows/local/star_align/main'
include { prepare_genome_workflow as PREPARE_GENOME } from './subworkflows/local/prepare_genome/main.nf'
include { initial_qc_workflow as INITIAL_QC }         from './subworkflows/local/initial_qc/main'
include { read_samples_workflow as READ_SAMPLES }     from './subworkflows/local/read_samples/main'





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
            fastqc_raw = Channel.empty()
            fastqvalidator = Channel.empty()
            cutadapt_reads = Channel.empty()
            cutadapt_log = Channel.empty()
            fastqc_trimmed = Channel.empty()
            bbtools_ihist = Channel.empty()
            fqscreen_1_txt = Channel.empty()
            fqscreen_1_png = Channel.empty()
            fqscreen_2_txt = Channel.empty()
            fqscreen_2_png = Channel.empty()

            star_pass1_sj = Channel.empty()
            star_pass1_log = Channel.empty()
            star_sjdb = Channel.empty()
            star_pass2_log = Channel.empty()
            star_pass2_sj = Channel.empty()
            star_pass2_reads_per_gene = Channel.empty()
            star_pass2_bam = Channel.empty()
            star_pass2_transcript_bam = Channel.empty()
        } else {
            log.info "Genome preparation complete. Continuing with workflow."
            prepare_genome_conf = Channel.empty() // dont save genome copy if not in build mode

            // Read samplesheet and emit channel of reads
            READ_SAMPLES(params.input)
            READ_SAMPLES.out.reads.set { ch_reads }

            // Initial QC, trimming, and FastQ Screen steps
            INITIAL_QC(ch_reads)
            fastqc_raw    = INITIAL_QC.out.fastqc_raw
            fastqvalidator = INITIAL_QC.out.fastqvalidator
            cutadapt_reads = INITIAL_QC.out.cutadapt_reads
            cutadapt_log   = INITIAL_QC.out.cutadapt_log
            fastqc_trimmed = INITIAL_QC.out.fastqc_trimmed
            bbtools_ihist  = INITIAL_QC.out.bbtools_ihist
            fqscreen_1_txt = INITIAL_QC.out.fqscreen_1_txt
            fqscreen_1_png = INITIAL_QC.out.fqscreen_1_png
            fqscreen_2_txt = INITIAL_QC.out.fqscreen_2_txt
            fqscreen_2_png = INITIAL_QC.out.fqscreen_2_png

            // STAR alignment steps ----------------------------------------------------------

            STAR_ALIGN(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf
            )
            star_pass1_sj            = STAR_ALIGN.out.pass1_sj
            star_pass1_log           = STAR_ALIGN.out.pass1_log
            star_sjdb                = STAR_ALIGN.out.sjdb
            star_pass2_log           = STAR_ALIGN.out.pass2_log
            star_pass2_sj            = STAR_ALIGN.out.pass2_sj
            star_pass2_reads_per_gene = STAR_ALIGN.out.pass2_reads_per_gene
            star_pass2_bam           = STAR_ALIGN.out.pass2_bam
            star_pass2_transcript_bam = STAR_ALIGN.out.pass2_transcript_bam
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
        prepare_genome_conf = prepare_genome_conf

        fastqc_raw    = fastqc_raw
        fastqvalidator = fastqvalidator
        cutadapt_reads = cutadapt_reads
        cutadapt_log   = cutadapt_log
        fastqc_trimmed = fastqc_trimmed
        bbtools_ihist  = bbtools_ihist
        fqscreen_1_txt = fqscreen_1_txt
        fqscreen_1_png = fqscreen_1_png
        fqscreen_2_txt = fqscreen_2_txt
        fqscreen_2_png = fqscreen_2_png

        star_pass1_sj            = star_pass1_sj
        star_pass1_log           = star_pass1_log
        star_sjdb                = star_sjdb
        star_pass2_log           = star_pass2_log
        star_pass2_sj            = star_pass2_sj
        star_pass2_reads_per_gene = star_pass2_reads_per_gene
        star_pass2_bam           = star_pass2_bam
        star_pass2_transcript_bam = star_pass2_transcript_bam
}

output {
    prepare_genome_conf { path { file -> "genome/" } }
    fastqc_raw { path { meta, file -> "fastqc/raw/" } }
    fastqvalidator { path { meta, file -> "fastqvalidator/${meta.id}/" } }
    cutadapt_reads { path { meta, reads -> "cutadapt/${meta.id}/" } }
    cutadapt_log { path { meta, log -> "cutadapt/${meta.id}/" } }
    fastqc_trimmed { path { meta, file -> "fastqc/trimmed/" } }
    bbtools_ihist { path { meta, ihist -> "bbtools/${meta.id}/" }}
    fqscreen_1_txt  { path { meta, file  -> "FQscreen/"                  } }
    fqscreen_1_png  { path { meta, file  -> "FQscreen/"                  } }
    fqscreen_2_txt  { path { meta, file  -> "FQscreen2/"                 } }
    fqscreen_2_png  { path { meta, file  -> "FQscreen2/"                 } }

    star_pass1_sj { path { meta, file -> 'STAR_files/pass1/' } }
    star_pass1_log { path { meta, file -> 'STAR_files/pass1/' } }
    star_sjdb { path { file -> 'STAR_files/pass1/' } }
    star_pass2_log { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_sj { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_reads_per_gene { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_bam { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_transcript_bam { path { meta, file -> 'bams/' } }
}
