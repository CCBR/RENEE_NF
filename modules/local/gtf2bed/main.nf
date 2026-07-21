process GTF2BED {
    tag { gtf }
    label 'process_single'

    container "${params.containers.base}"

    input:
        tuple val(meta), path(gtf)

    output:
        path("${gtf.baseName}.bed"), emit: bed

    script:
    """
    # remove comment lines and convert GTF to BED format.
    # looks for transcript_id for test case
    grep -v "^#" ${gtf} | grep 'transcript_id' | gtf2bed > ${gtf.baseName}.bed

    """
    stub:
    """
    touch ${gtf.baseName}.bed
    """
}
