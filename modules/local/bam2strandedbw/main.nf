process BAM2STRANDEDBW {
    tag { meta.id }
    label 'process_high'

    container 'nciccbr/ccbr_bam2strandedbw:v0.0.1'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.fwd.bw"), emit: fwd_bw
        tuple val(meta), path("*.rev.bw"), emit: rev_bw
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed "s/samtools //g"'), emit: versions_samtools, topic: versions
        tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed "s/bedtools //g"'), emit: versions_bedtools, topic: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix     ?: "${meta.id}"
    def swap       = task.ext.swap_strands ?: false

    # TODO - Add support for strandedness detection (e.g. using RSeQC infer_experiment.py)
    if (meta.single_end) {
        """
        # Extract chromosome sizes from BAM header
        samtools view -H ${bam} \\
            | grep "^@SQ" \\
            | cut -f2,3 \\
            | sed 's/SN://g' \\
            | sed 's/LN://g' \\
            > ${bam}.genome

        # Paired flags for SE:
        #   forward strand: reads on reverse complement (-f 16 = read mapped on reverse strand)
        #   reverse strand: reads on forward strand (-F 16 = read NOT on reverse strand)
        samtools view -b -f 16 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.fwd.bg
        bedSort ${prefix}.fwd.bg ${prefix}.fwd.bg

        samtools view -b -F 16 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.rev.bg
        bedSort ${prefix}.rev.bg ${prefix}.rev.bg

        bedGraphToBigWig ${prefix}.fwd.bg ${bam}.genome ${prefix}.fwd.bw
        bedGraphToBigWig ${prefix}.rev.bg ${bam}.genome ${prefix}.rev.bw

        rm -f ${prefix}.fwd.bg ${prefix}.rev.bg ${bam}.genome

        # Swap fwd/rev if requested (non-dUTP / FIRST_READ_TRANSCRIPTION_STRAND)
        if [ "${swap}" = "true" ]; then
            mv ${prefix}.fwd.bw ${prefix}.fwd.bw.tmp
            mv ${prefix}.rev.bw ${prefix}.fwd.bw
            mv ${prefix}.fwd.bw.tmp ${prefix}.rev.bw
        fi
        """
    } else {
        """
        # Extract chromosome sizes from BAM header
        samtools view -H ${bam} \\
            | grep "^@SQ" \\
            | cut -f2,3 \\
            | sed 's/SN://g' \\
            | sed 's/LN://g' \\
            > ${bam}.genome

        # Compute forward strand bedgraphs in parallel (dUTP convention):
        #   fwd1: read2, forward  (-f 128 -F 16)
        #   fwd2: read1, reverse  (-f 80)
        samtools view -b -f 128 -F 16 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.fwd1.bg &
        samtools view -b -f 80 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.fwd2.bg &

        # Compute reverse strand bedgraphs in parallel:
        #   rev1: read2, reverse  (-f 144)
        #   rev2: read1, forward  (-f 64 -F 16)
        samtools view -b -f 144 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.rev1.bg &
        samtools view -b -f 64 -F 16 ${bam} \\
            | bedtools genomecov -bg -split -ibam stdin -g ${bam}.genome \\
            > ${prefix}.rev2.bg &
        wait

        bedSort ${prefix}.fwd1.bg ${prefix}.fwd1.bg
        bedSort ${prefix}.fwd2.bg ${prefix}.fwd2.bg
        bedSort ${prefix}.rev1.bg ${prefix}.rev1.bg
        bedSort ${prefix}.rev2.bg ${prefix}.rev2.bg

        # Merge sub-strands into single fwd/rev bedgraphs
        bedtools unionbedg -i ${prefix}.fwd1.bg ${prefix}.fwd2.bg \\
            | awk -F'\\t' -v OFS='\\t' '{sum=0; for(i=4;i<=NF;i++) sum+=\$i; print \$1,\$2,\$3,sum}' \\
            > ${prefix}.fwd.bg
        bedtools unionbedg -i ${prefix}.rev1.bg ${prefix}.rev2.bg \\
            | awk -F'\\t' -v OFS='\\t' '{sum=0; for(i=4;i<=NF;i++) sum+=\$i; print \$1,\$2,\$3,sum}' \\
            > ${prefix}.rev.bg

        bedSort ${prefix}.fwd.bg ${prefix}.fwd.bg
        bedSort ${prefix}.rev.bg ${prefix}.rev.bg

        bedGraphToBigWig ${prefix}.fwd.bg ${bam}.genome ${prefix}.fwd.bw
        bedGraphToBigWig ${prefix}.rev.bg ${bam}.genome ${prefix}.rev.bw

        rm -f ${prefix}.fwd*.bg ${prefix}.rev*.bg ${bam}.genome

        # Swap fwd/rev if requested (non-dUTP / FIRST_READ_TRANSCRIPTION_STRAND)
        if [ "${swap}" = "true" ]; then
            mv ${prefix}.fwd.bw ${prefix}.fwd.bw.tmp
            mv ${prefix}.rev.bw ${prefix}.fwd.bw
            mv ${prefix}.fwd.bw.tmp ${prefix}.rev.bw
        fi
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fwd.bw
    touch ${prefix}.rev.bw
    """
}
