nextflow.enable.dsl = 2

// Plugins
// Modules
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include {FASTQC as FASTQC_RAW}                   from './modules/local/fastqc'
include {FASTQC as FASTQC_TRIMMED}               from './modules/local/fastqc'
include {BBTOOLS_BBMERGE}                        from './modules/local/bbtools'
include {CUTADAPT}                               from './modules/CCBR/cutadapt'
include {FASTQSCREEN_FASTQSCREEN as FASTQ_SCREEN_1} from './modules/nf-core/fastqscreen/fastqscreen/main.nf'
include {FASTQSCREEN_FASTQSCREEN as FASTQ_SCREEN_2} from './modules/nf-core/fastqscreen/fastqscreen/main.nf'

// Subworkflows
include {star_align as STAR_ALIGN}               from './subworkflows/local/star_align/main'
include {validate_fastqs as VALIDATE_FASTQS}     from './subworkflows/local/validate_fastqs/main'
include { prepare_genome as PREPARE_GENOME }     from './subworkflows/local/prepare_genome/main.nf'





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

workflow MAKE_REFERENCE {
    PREPARE_GENOME()
}

workflow {
    main:
        LOG()
        validateParameters()

        ch_sjdb_placeholder = Channel.value(file(params.sjdb_placeholder_tab, checkIfExists: true))

        ch_reads = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def has_fastq_2 = row.fastq_2 && row.fastq_2.toString().trim()
                def meta = [
                    id       : "${row.sample}_${row.replicate}",
                    sample   : row.sample,
                    replicate: row.replicate,
                    single_end   : !has_fastq_2
                ]

                def reads = [file(row.fastq_1)]
                if (has_fastq_2) {
                    reads << file(row.fastq_2)
                }

                tuple(meta, reads)
            }
        // Split each sample read list into one fastq per emitted tuple for validation.
        individual_fastq_ch = ch_reads.transpose()


        // Sample validation gate
        VALIDATE_FASTQS(individual_fastq_ch, ch_reads)
        ch_validated_reads = VALIDATE_FASTQS.out.reads

        // ch_validated_reads.view()

        // QC and trimming steps
        FASTQC_RAW(ch_validated_reads)

        CUTADAPT(ch_validated_reads)

        FASTQC_TRIMMED(CUTADAPT.out.reads)
        BBTOOLS_BBMERGE(CUTADAPT.out.reads)

        ch_fqscreen_1_txt = Channel.empty()
        ch_fqscreen_1_png = Channel.empty()
        ch_fqscreen_2_txt = Channel.empty()
        ch_fqscreen_2_png = Channel.empty()


        // FastQ Screen steps
        if (!params.fastq_screen_db_dir) {
            log.warn "No FastQ Screendatabase directory provided. FastQ Screen will be skipped."
        } else {
            ch_fqscreen_db_dir = Channel.value(file(params.fastq_screen_db_dir))

            if (params.fastq_screen_conf) {
                FASTQ_SCREEN_1(CUTADAPT.out.reads, file(params.fastq_screen_conf), ch_fqscreen_db_dir)
                ch_fqscreen_1_txt = FASTQ_SCREEN_1.out.txt
                ch_fqscreen_1_png = FASTQ_SCREEN_1.out.png
            }

            if (params.fastq_screen_conf2) {
                FASTQ_SCREEN_2(CUTADAPT.out.reads, file(params.fastq_screen_conf2), ch_fqscreen_db_dir)
                ch_fqscreen_2_txt = FASTQ_SCREEN_2.out.txt
                ch_fqscreen_2_png = FASTQ_SCREEN_2.out.png
            }
        }


        // STAR alignment steps ----------------------------------------------------------

        PREPARE_GENOME()


        STAR_ALIGN(
            CUTADAPT.out.reads,
            ch_sjdb_placeholder,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.genes_gtf
            )

        workflow.onComplete = {
            if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
                def message = Utils.spooker(workflow)
                if (message) {
                    println message
                }
            }
        }

    publish:
        fastqc_raw = FASTQC_RAW.out.html.mix(FASTQC_RAW.out.zip)
        fastqvalidator = VALIDATE_FASTQS.out.logs
        cutadapt_reads = CUTADAPT.out.reads
        cutadapt_log = CUTADAPT.out.log
        fastqc_trimmed = FASTQC_TRIMMED.out.html.mix(FASTQC_TRIMMED.out.zip)
        bbtools_ihist = BBTOOLS_BBMERGE.out.ihist
        fqscreen_1_txt = ch_fqscreen_1_txt
        fqscreen_1_png = ch_fqscreen_1_png
        fqscreen_2_txt = ch_fqscreen_2_txt
        fqscreen_2_png = ch_fqscreen_2_png

        star_pass1_sj = STAR_ALIGN.out.pass1_sj
        star_pass1_log = STAR_ALIGN.out.pass1_log
        star_sjdb = STAR_ALIGN.out.sjdb
        star_pass2_log = STAR_ALIGN.out.pass2_log
        star_pass2_sj = STAR_ALIGN.out.pass2_sj
        star_pass2_reads_per_gene = STAR_ALIGN.out.pass2_reads_per_gene
        star_pass2_bam = STAR_ALIGN.out.pass2_bam
        star_pass2_transcript_bam = STAR_ALIGN.out.pass2_transcript_bam
}

output {
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
