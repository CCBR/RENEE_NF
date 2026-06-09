include { FASTQVALIDATOR } from "../../../modules/local/fastqvalidator"

workflow validate_fastqs {
    take:
        ch_individual_fastqs
        ch_reads

    main:
        FASTQVALIDATOR(ch_individual_fastqs)

        ch_validator_check = FASTQVALIDATOR.out.result
            .map { meta, fastq, log, exitcode_file ->
                [
                    id      : meta.id,
                    fastq   : fastq,
                    exitcode: exitcode_file.text.trim() as Integer,
                    log     : log
                ]
            }
            .collect()
            .map { results ->
                def failed = results.findAll { entry -> entry.exitcode != 0 }

                if (failed) {
                    def failed_summary = failed.collect { entry ->
                        "${entry.id}\t${entry.fastq.baseName}"
                    }.join('\n')
                    error "FASTQVALIDATOR failed for FASTQ file(s):\n${failed_summary}"
                }

                true
            }

        // Hold downstream reads until all validator tasks are complete and checked.
        ch_validated_reads = ch_reads
            .combine(ch_validator_check)
            .map { meta, reads, validator -> tuple(meta, reads) }

        // make a channel of all the logs for downstream aggregation
        FASTQVALIDATOR.out.result
            .map { meta, fastq, log, exitcode_file -> tuple(meta, log) }
            .set { ch_fastqvalidator_logs }

    emit:
        reads = ch_validated_reads
        logs = ch_fastqvalidator_logs
}
