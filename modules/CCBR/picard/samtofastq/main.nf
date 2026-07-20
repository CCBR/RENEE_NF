process PICARD_SAMTOFASTQ {
    tag { meta.id }
    label 'process_medium'

    container 'nciccbr/ccbr_picard_2.27.5:v1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"),         emit: reads
    tuple val(meta), path("*_?.fastq.gz"),       emit: paired, optional: true
    tuple val(meta), path("*unpaired.fastq.gz"), emit: unpaired, optional: true
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = meta.single_end ? "--FASTQ ${prefix}_1.fastq" : "--FASTQ ${prefix}_1.fastq --SECOND_END_FASTQ ${prefix}_2.fastq --UNPAIRED_FASTQ ${prefix}.unpaired.fastq"

    if (!task.memory) {
        log.warn '[Picard SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072

    """
    picard \\
        -Xmx${avail_mem}M \\
        SamToFastq \\
        ${args} \\
        --INPUT ${bam} \\
        ${output}

    pigz -p ${task.cpus} *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SamToFastq --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1.fastq.gz versions.yml
    """
}
