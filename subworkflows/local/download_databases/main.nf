include { UNTAR as UNTAR_FASTQ_DB  } from '../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_KRAKEN_DB } from '../../../modules/nf-core/untar/main.nf'


workflow DOWNLOAD_DATABASES {

    main:

        // Kraken and FastQ have different file structures, so we need to untar them separately to format results correctly
        ch_all = Channel
            .fromList(params.database_urls.collect { db_name, url -> tuple(db_name, url) })

        ch_fastq_archives = ch_all
            .filter { db_name, url -> db_name != '20180907_standard_kraken2' }
            .map { db_name, url -> tuple([ id: db_name ], file(url)) }

        ch_kraken_archives = ch_all
            .filter { db_name, url -> db_name == '20180907_standard_kraken2' }
            .map { db_name, url -> tuple([ id: db_name ], file(url)) }

        UNTAR_FASTQ_DB(ch_fastq_archives)
        UNTAR_KRAKEN_DB(ch_kraken_archives)

    emit:
        fastq_screen_databases = UNTAR_FASTQ_DB.out.untar
        kraken_databases       = UNTAR_KRAKEN_DB.out.untar
}
