process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "nciccbr/ccbr_arriba_2.0.0:v0.0.1"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val star_use_sjdb
    path sjdb

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed "s/STAR_//"'), emit: versions_star, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('gawk'), eval("gawk --version | sed -n '1s/GNU Awk \\([0-9.]*\\).*/\\1/p'"), emit: versions_gawk, topic: versions

    tuple val(meta), path('*d.out.bam')                              , optional:true, emit: bam
    tuple val(meta), path("${prefix}.sortedByCoord.out.bam")         , optional:true, emit: bam_sorted
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam") , optional:true, emit: bam_sorted_aligned
    tuple val(meta), path('*toTranscriptome.out.bam')                , optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam')                 , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')                               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')                            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')                  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')                          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')                               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                                    , optional:true, emit: bedgraph
    tuple val(meta), path('*._STARgenome')                           , optional:true, emit: stargenome

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each{ read -> reads1 << read} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def sjdb_arg        = star_use_sjdb ? "--sjdbFileChrStartEnd $sjdb" : ''
    attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' 'SM:$prefix'"
    def out_sam_type    = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        --sjdbGTFfile $gtf \\
        $sjdb_arg \\
        $attrRG \\
        $args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    # is not always cleaned up by STAR, especially in case of errors, and can take up a lot of space
    rm -rf ${prefix}._STARtmp
    rm -rf ${prefix}._STARgenome
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg
    mkdir ${prefix}._STARgenome
    touch ${prefix}._STARgenome/exonGeTrInfo.tab ${prefix}._STARgenome/exonInfo.tab ${prefix}._STARgenome/geneInfo.tab ${prefix}._STARgenome/sjdbInfo.txt ${prefix}._STARgenome/sjdbList.fromGTF.out.tab ${prefix}._STARgenome/sjdbList.out.tab ${prefix}._STARgenome/transcriptInfo.tab
    """
}
