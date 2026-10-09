workflow CHECK_INPUT {
    take:
        samplesheet  // path: CSV samplesheet

    main:
        ch_reads = Channel
            .fromPath(samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def has_fastq_2 = row.fastq_2 && row.fastq_2.toString().trim()
                def meta = [
                    id       : "${row.sample}_${row.replicate}",
                    sample   : row.sample,
                    replicate: row.replicate,
                    single_end   : !has_fastq_2
                ]

                def projectDirPath = { path ->
                    file(path.toString().replace('${projectDir}', projectDir.toString()))
                }
                def reads = [projectDirPath(row.fastq_1)]
                if (has_fastq_2) {
                    reads << projectDirPath(row.fastq_2)
                }

                if (params.small_rna && has_fastq_2) {
                    error "The small_rna option is only supported with single-end data, but sample ${meta.id} is paired-end."
                }

                tuple(meta, reads)
            }

    emit:
        reads = ch_reads  // channel: [ meta, [ reads ] ]
}
