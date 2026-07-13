// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process FUSIONARRIBA {
    tag "$meta.id"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "nciccbr/ccbr_arriba_2.0.0:v0.0.1"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(index)
    path(blacklist)

    output:
    tuple val(meta), path("*.fusions.tsv")          , emit: fusions
    tuple val(meta), path("*.fusions.discarded.tsv"), emit: fusions_fail
    tuple val(meta), path("*.Aligned.out.bam")       , emit: bam
    tuple val(meta), path("*Log.final.out")          , emit: log_final
    tuple val(meta), path("*Log.out")                , emit: log_out
    tuple val(meta), path("*Log.progress.out")       , emit: log_progress
    tuple val("${task.process}"), val('star'),   eval('STAR --version | sed "s/STAR_//"'),                                   emit: versions_star,   topic: versions
    tuple val("${task.process}"), val('arriba'), eval('arriba -h | grep \'Version:\' 2>&1 |  sed \'s/Version:\\s//\''), emit: versions_arriba, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        error "FUSIONARRIBA requires paired-end reads, but single-end data was provided for sample: ${meta.id}"
    }

    """
    # Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
    readlength=\$(
        zcat ${reads[0]} | \\
        awk -v maxlen=100 'NR%4==2 {if (length(\$1) > maxlen+0) maxlen=length(\$1)}; \\
        END {print maxlen-1}'
    )

    STAR --runThreadN ${task.cpus} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang \${readlength} \\
        --genomeDir ${index} \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outStd BAM_Unsorted \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --outFilterMultimapNmax 50 \\
        --peOverlapNbasesMin 10 \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 \\
        --chimOutType WithinBAM HardClip \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 \\
        --chimMultimapNmax 50 \\
        --twopassMode Basic \\
        --outFileNamePrefix ${prefix}. \\
        ${args} \\
    | tee ${prefix}.Aligned.out.bam | \\
    arriba -x /dev/stdin \\
        -o ${prefix}.fusions.tsv \\
        -O ${prefix}.fusions.discarded.tsv \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b ${blacklist} \\
        ${args2}

    # Temporary STAR directories
    # are not always cleaned up by STAR, especially in case of errors, and can take up a lot of space
    rm -rf ${prefix}._STARtmp
    rm -rf ${prefix}._STARgenome
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        error "FUSIONARRIBA requires paired-end reads, but single-end data was provided for sample: ${meta.id}"
    }

    """
    echo $args
    touch ${prefix}.fusions.tsv
    touch ${prefix}.fusions.discarded.tsv
    touch ${prefix}.Aligned.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    """
}
