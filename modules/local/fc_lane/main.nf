process FC_LANE {
    tag { meta.id }
    label 'process_single'
    label 'qc'

    conda "${moduleDir}/environment.yml"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.fastq.info.txt"), emit: fqinfo
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | cut -d " " -f 2'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_flowcell_lanes.py \
        ${fastq} \
        ${meta.id} \
        > ${prefix}.fastq.info.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample_name\\ttotal_read_pairs\\tflowcell_ids\\tlanes\\tflowcell_lanes\\tmd5_checksum\\n' \
        > ${prefix}.fastq.info.txt
    printf '${meta.id}\\t0\\tNA\\tNA\\tNA\\tNA\\n' \
        >> ${prefix}.fastq.info.txt
    """
}
