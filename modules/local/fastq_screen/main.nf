process FASTQ_SCREEN {
    tag { meta.id }
    label 'process_high'

    container "nciccbr/ccbr_fastq_screen_0.13.0:v2.0"

    input:
        tuple val(meta), path(reads)
        path fastq_screen_config

    output:
        tuple val(meta), path("*_screen.txt"), emit: txt
        tuple val(meta), path("*_screen.png"), emit: png
        tuple val("${task.process}"), val('fastq_screen'), eval('fastq_screen --version 2>&1 | head -1'), emit: versions_fastq_screen, topic: versions

    script:
    /*
     * Quality-control step to screen for different sources of contamination.
     * FastQ Screen compares sequencing data to a set of reference genomes to
     * determine if there is contamination. Call this process twice (aliased) in
     * the workflow to reproduce the dual-screen behaviour (multi-organism and
     * vector/rRNA) that the Snakemake pipeline performed in a single rule.
     *
     * @Input:  Trimmed FastQ file(s) — single-end (R1 only) or paired-end (R1 + R2)
     * @Output: FastQ Screen .txt and .png reports
     */
    def reads_arg = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    fastq_screen \\
        --conf ${fastq_screen_config} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        ${reads_arg}
    """

    stub:
    def r1_stem = reads[0].name
        .replaceAll(/\.fastq\.gz$/, '')
        .replaceAll(/\.fastq$/, '')
    def r2_stem = meta.single_end ? null
        : reads[1].name
            .replaceAll(/\.fastq\.gz$/, '')
            .replaceAll(/\.fastq$/, '')
    """
    touch ${r1_stem}_screen.txt
    touch ${r1_stem}_screen.png
    ${r2_stem ? "touch ${r2_stem}_screen.txt" : ''}
    ${r2_stem ? "touch ${r2_stem}_screen.png" : ''}
    """
}
