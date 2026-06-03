
process FASTQC {
    tag "${meta.id}:${meta.qc_stage ?: 'raw'}"

    container "nciccbr/ccrgb_qctools:v4.0"
    publishDir { [ path: "${params.outputDir}/fastqc/${meta.qc_stage ?: 'raw'}", mode: params.publish_dir_mode ] }
    input:
        tuple val(meta), path(fastq)
    output:
        tuple val(meta), path("*fastqc.html"), emit: html
        tuple val(meta), path("*fastqc.zip"),  emit: zip
        tuple val("${task.process}"), val('fastqc'), eval('fastqc --version | sed "/FastQC v/!d; s/.*v//"'), emit: versions_fastqc, topic: versions

    script:
    """
    fastqc \
        $fastq \
        -t $task.cpus \
        -o .
    """

    stub:
    """
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip
    """
}
