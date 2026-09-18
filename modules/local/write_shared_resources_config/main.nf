process WRITE_SHARED_RESOURCES_CONFIG {
    label 'process_single'
    container "${params.containers.base}"

    input:
        path(kraken_db_dir)
        path(fastq_screen_conf1)
        path(fastq_screen_conf2)

    output:
        path("shared_resources.config"), emit: conf

    script:
    // Sibling of shared_resources/ itself, mirroring WRITE_GENOME_CONFIG's
    // <genome name>.config sitting alongside genome/<genome name>/. Absolute,
    // so a later analysis run's `-c` keeps working regardless of its own
    // --outputDir -- matches how the databases this points at are actually
    // published (see Utils.sharedResourcesDir / main.nf's output {} block).
    def shared_resources_dir = Utils.sharedResourcesDir(params)
    """
    {
        echo 'params {'
        echo '    kraken2_db_dir     = "${shared_resources_dir}/${kraken_db_dir.name}"'
        echo '    fastq_screen_conf  = "${shared_resources_dir}/fastq_screen_db/${fastq_screen_conf1.name}"'
        echo '    fastq_screen_conf2 = "${shared_resources_dir}/fastq_screen_db/${fastq_screen_conf2.name}"'
        echo '}'
    } > shared_resources.config
    """

    // No stub block: this just writes 3 lines of text from filenames already
    // known to Nextflow -- cheap enough to always run for real, even under
    // -stub-run, unlike the heavier downloads/untars it runs after.
}
