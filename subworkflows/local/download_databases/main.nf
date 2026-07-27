include { UNTAR as UNTAR_DATABASE } from '../../../modules/nf-core/untar/main.nf'


workflow DOWNLOAD_DATABASES {

    main:

        /*
         * Convert each configured database URL into the nf-core UNTAR input shape:
         * tuple val(meta), path(archive).
         */
        ch_database_archives = Channel
            .fromList(params.database_urls.collect { db_name, url -> tuple(db_name, url) })
            .map { db_name, url ->
                def meta = [
                    id      : db_name,
                    db_group: db_name == '20180907_standard_kraken2' ? null : 'fastq_screen_db'
                ]

                tuple(meta, file(url))
            }

        UNTAR_DATABASE(ch_database_archives)

    emit:
        databases = UNTAR_DATABASE.out.untar
}
