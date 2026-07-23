include { FASTQC as FASTQC_RAW }                        from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_TRIMMED }                    from '../../../modules/local/fastqc'
include { BBTOOLS_BBMERGE }                             from '../../../modules/local/bbtools'
include { CUTADAPT }                                    from '../../../modules/CCBR/cutadapt'
include { FASTQSCREEN_FASTQSCREEN as FASTQ_SCREEN_1 }   from '../../../modules/nf-core/fastqscreen/fastqscreen/main.nf'
include { FASTQSCREEN_FASTQSCREEN as FASTQ_SCREEN_2 }   from '../../../modules/nf-core/fastqscreen/fastqscreen/main.nf'
include { VALIDATE_FASTQS }                             from '../validate_fastqs/main'
include { KRAKEN2_KRAKEN2 }                             from '../../../modules/nf-core/kraken2/kraken2/main'
include { UNTAR as KRAKEN2_DB_DOWNLOAD }                                       from '../../../modules/nf-core/untar'

workflow INITIAL_QC {
    take:
        ch_reads  // channel: [ meta, [ reads ] ]

    main:
        // Split each sample read list into one fastq per emitted tuple for validation.
        individual_fastq_ch = ch_reads.transpose()

        // Sample validation gate
        VALIDATE_FASTQS(individual_fastq_ch, ch_reads)
        ch_validated_reads = VALIDATE_FASTQS.out.reads

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
            log.warn "No FastQ Screen database directory provided. FastQ Screen will be skipped."
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

        // Kraken2 taxonomic classification step
        ch_kraken2_db_dir = Channel.empty()

        // First check if db is present
        // Then if a link to download the db is provided
        // finally if no db is provided, log a warning and skip the step
        if (params.kraken2_db_dir){
            ch_kraken2_db_dir = Channel.value(file(params.kraken2_db_dir))

        } else if (params.kraken2_db_url) {
            // Download the Kraken2 database from the provided URL
            // and unzip it
            ch_kraken2_url = Channel.of([[id:'kraken2_db'],  file(params.kraken2_db_url)])
            KRAKEN2_DB_DOWNLOAD(ch_kraken2_url)
            ch_kraken2_db_dir = KRAKEN2_DB_DOWNLOAD.out.untar
                .map { meta, path -> path }
                .first()

        } else {
            log.warn "No Kraken2 database directory, or download URL provided. Kraken2 will be skipped."
        }
        // this should skip with an empty channel if no db is provided
        // TODO: test that
        KRAKEN2_KRAKEN2(ch_validated_reads, ch_kraken2_db_dir, params.kraken_save_output_fastqs, params.kraken_save_reads_assignment)


    emit:
        trimmed_reads   = CUTADAPT.out.reads

        fastqvalidator  = VALIDATE_FASTQS.out.logs

        fastqc_raw      = FASTQC_RAW.out.html.mix(FASTQC_RAW.out.zip)

        cutadapt_reads  = CUTADAPT.out.reads
        cutadapt_log    = CUTADAPT.out.log

        fastqc_trimmed  = FASTQC_TRIMMED.out.html.mix(FASTQC_TRIMMED.out.zip)

        bbtools_ihist   = BBTOOLS_BBMERGE.out.ihist

        fqscreen_1_txt  = ch_fqscreen_1_txt
        fqscreen_1_png  = ch_fqscreen_1_png
        fqscreen_2_txt  = ch_fqscreen_2_txt
        fqscreen_2_png  = ch_fqscreen_2_png

        kraken2_db_dir  = ch_kraken2_db_dir

        kraken2_report  = KRAKEN2_KRAKEN2.out.report
        kraken2_classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
        kraken2_krona_html = KRAKEN2_KRAKEN2.out.krona_html
}
