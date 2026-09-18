process CREATE_FASTQ_SCREEN_CONF {
    tag "fastq_screen_conf"
    label 'process_single'
    container "${params.containers.base}"

    input:
        path(conf1_template)
        path(conf2_template)

    output:
        path("fastq_screen_p1.conf"), emit: conf1
        path("fastq_screen_p2.conf"), emit: conf2

    when:
        task.ext.when == null || task.ext.when

    script:
    // Mirrors RENEE/workflow/rules/build.smk rule fqscreen_conf: repoint the
    // DATABASE paths of the checked-in FastQ Screen conf templates
    // (assets/fastq_screen_p{1,2}.conf.template) at the shared resources
    // directory this build publishes to (always <outputDir>/shared_resources,
    // see Utils.sharedResourcesDir). The templates are local copies of
    // https://hpc.nih.gov/~OpenOmics/common/fastq_screen_p{1,2}.conf with
    // their DATABASE paths already pointed at a prior local build; that
    // hard-coded path is reused below as the placeholder to substitute, so
    // this no longer needs outbound network access from compute nodes.
    def db_root = Utils.sharedResourcesDir(params)
    def placeholder = '/projectnb/wax-es/alecs/renee/RENEE_NF/build_results2/shared_resources'
    """
    sed 's@${placeholder}@${db_root}@g' ${conf1_template} > fastq_screen_p1.conf
    sed 's@${placeholder}@${db_root}@g' ${conf2_template} > fastq_screen_p2.conf
    """

    stub:
    """
    touch fastq_screen_p1.conf fastq_screen_p2.conf
    """
}
