process BBTOOLS_BBMERGE {
    tag { meta.id }
    label 'process_high'

    container 'nciccbr/ccbr_bbtools_38.87:v0.0.1'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*_insert_sizes.txt'), emit: ihist
        tuple val("${task.process}"), val('bbtools'), eval('bbmerge.sh --version 2>&1 | head -n 1'), emit: versions_bbtools, topic: versions

    when:
        (task.ext.when == null || task.ext.when) && meta.layout == 'paired'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def qin = task.ext.qin ?: 33
    def read1 = reads[0]
    def read2 = reads[1]
    def xmxGb = task.memory ? Math.max(1, task.memory.toGiga() as Integer) : 64

    """
    bbtools bbmerge-auto \
        in1=${read1} \
        in2=${read2} \
        qin=${qin} \
        ihist=${prefix}_insert_sizes.txt \
        k=62 \
        extend2=200 \
        rem \
        ecct \
        -Xmx${xmxGb}g
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_insert_sizes.txt
    """
}
