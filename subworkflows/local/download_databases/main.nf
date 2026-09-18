include { UNTAR as UNTAR_FASTQ_DB      } from '../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_KRAKEN_DB     } from '../../../modules/nf-core/untar/main.nf'
include { DOWNLOAD_FASTQ_SCREEN_CONF   } from '../../../modules/local/download_fastq_screen_conf/main.nf'
include { ARRIBA_DOWNLOAD              } from '../../../modules/nf-core/arriba/download/main.nf'


workflow DOWNLOAD_DATABASES {

    take:
        shared_resources_dir

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

        // FastQ Screen conf files (not archived with the DBs above); mirrors
        // RENEE/workflow/rules/build.smk rule fqscreen_conf. Templates are
        // checked-in local assets rather than a runtime download -- see
        // modules/local/download_fastq_screen_conf/main.nf.
        DOWNLOAD_FASTQ_SCREEN_CONF(
            shared_resources_dir,
            file("${projectDir}/assets/fastq_screen_p1.conf.template"),
            file("${projectDir}/assets/fastq_screen_p2.conf.template")
        )

        // Arriba fusion-calling reference database. RENEE (classic) never
        // downloads this -- it's provisioned here as a genome-agnostic
        // shared resource instead (one tarball covers every genome build).
        // Passing an empty genome value keeps ARRIBA_DOWNLOAD's glob outputs
        // unfiltered, i.e. every genome's files, not just one.
        ARRIBA_DOWNLOAD('')

    emit:
        fastq_screen_databases = UNTAR_FASTQ_DB.out.untar
        kraken_databases       = UNTAR_KRAKEN_DB.out.untar
        fastq_screen_conf1     = DOWNLOAD_FASTQ_SCREEN_CONF.out.conf1
        fastq_screen_conf2     = DOWNLOAD_FASTQ_SCREEN_CONF.out.conf2
        arriba_database        = ARRIBA_DOWNLOAD.out.blacklist
            .mix(ARRIBA_DOWNLOAD.out.cytobands)
            .mix(ARRIBA_DOWNLOAD.out.protein_domains)
            .mix(ARRIBA_DOWNLOAD.out.known_fusions)
}
