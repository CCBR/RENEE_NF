process BUILD_RSEM_REF {
    label 'process_high'
    container "nciccbr/ccbr_rsem_1.3.3:v1.0"

    input:
        path(fasta)
        path(gtf)
        val(genome_name)

    output:
        path("rsemref"), emit: rsem_ref

    when:
        task.ext.when == null || task.ext.when

    script:
    // Ports RENEE (classic) workflow/rules/build.smk rule `rsem`. Prefix
    // (rsemref/<genome_name>) is written into the generated genome config's
    // rsem_ref field by WRITE_GENOME_CONFIG -- rsem-calculate-expression
    // takes it as a --reference prefix, not a single file.
    """
    mkdir -p rsemref
    rsem-prepare-reference -p ${task.cpus} --gtf ${gtf} ${fasta} rsemref/${genome_name}
    rsem-generate-ngvector rsemref/${genome_name}.transcripts.fa rsemref/${genome_name}.transcripts
    """

    stub:
    """
    mkdir -p rsemref
    touch rsemref/${genome_name}.grp rsemref/${genome_name}.transcripts.fa
    """
}
