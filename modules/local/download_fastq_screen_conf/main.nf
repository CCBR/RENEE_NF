process DOWNLOAD_FASTQ_SCREEN_CONF {
    tag "fastq_screen_conf"
    label 'process_single'
    container "${params.containers.base}"

    input:
        val(shared_resources_dir)

    output:
        path("fastq_screen_p1.conf"), emit: conf1
        path("fastq_screen_p2.conf"), emit: conf2

    when:
        task.ext.when == null || task.ext.when

    script:
    // Mirrors RENEE/workflow/rules/build.smk rule fqscreen_conf: download the
    // FastQ Screen conf files and repoint their DATABASE paths at the shared
    // resources directory the databases are published to.
    def db_root = "${shared_resources_dir}".replaceAll('/+$', '')
    """
    wget https://hpc.nih.gov/~OpenOmics/common/fastq_screen_p1.conf -O fastq_screen_p1.conf
    wget https://hpc.nih.gov/~OpenOmics/common/fastq_screen_p2.conf -O fastq_screen_p2.conf
    sed -i 's@/data/OpenOmics/references/common@${db_root}@g' fastq_screen_p1.conf
    sed -i 's@/data/OpenOmics/references/common@${db_root}@g' fastq_screen_p2.conf
    """

    stub:
    """
    touch fastq_screen_p1.conf fastq_screen_p2.conf
    """
}
